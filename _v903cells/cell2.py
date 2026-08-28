def random_unit_vector(n=1, rng=rng):
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    return v[0] if n == 1 else v


def wrap_deg(x):
    return (x + 180.0) % 360.0 - 180.0


class CollectiveNanocubeMC:
    """Rigid nanocube cluster with independent positions, bodies, and dipoles."""

    def __init__(self, L_nm, phys=PHYS_V602, anis_model="cubic_first_raw"):
        self.L_nm = float(L_nm)
        self.L = self.L_nm * 1e-9
        self.vol = self.L**3
        self.phys = dict(phys)
        self.Ms = phys["Ms"]
        self.mu_mag = self.Ms * self.vol
        self.Bvec = np.array([phys["B_field"], 0.0, 0.0])
        self.Bhat = self.Bvec / np.linalg.norm(self.Bvec)
        self.K = phys["K_ani"]
        self.Aham = phys["Hamaker"]
        self.kT = const.k * T_K
        self.Cmag = const.mu_0 / (4.0 * np.pi)
        self.gap_m = phys["gap_nm"] * 1e-9
        self.gap_nm = phys["gap_nm"]
        self.k_stiff = phys["k_stiff"]
        self.r_round = min(phys["roundness_nm"] * 1e-9, self.L / 2.0 - 1e-12)
        self.anis_model = anis_model

        Nv = phys["N_voxel"]
        step = self.L / Nv
        lin = np.linspace(-self.L / 2.0 + step / 2.0, self.L / 2.0 - step / 2.0, Nv)
        gx, gy, gz = np.meshgrid(lin, lin, lin)
        self.vox_ref = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])
        self.vox_vol = step**3
        self.scv = -(self.Aham / np.pi**2) * self.vox_vol**2

        self.e111 = np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0)
        self.body_u = np.array([1.0, -1.0, 0.0]) / np.sqrt(2.0)
        self.field_hat = np.array([1.0, 0.0, 0.0])
        self.Rbase = R.align_vectors([[1.0, 0.0, 0.0]], [self.e111])[0].as_matrix()

    def rhombo_basis(self, a_nm, alpha_deg):
        a = a_nm * 1e-9
        al = np.radians(alpha_deg)
        cp = np.sqrt(max((2.0 * np.cos(al) + 1.0) / 3.0, 0.0))
        sp = np.sqrt(max(1.0 - cp**2, 0.0))
        return np.array(
            [
                [cp, sp, 0.0],
                [cp, -0.5 * sp, 0.866 * sp],
                [cp, -0.5 * sp, -0.866 * sp],
            ]
        ) * a

    def rotation_parts(self, beta_deg):
        Ry = R.from_euler("y", beta_deg, degrees=True).as_matrix()
        chain_hat = Ry @ self.field_hat
        return Ry, chain_hat

    def lattice_basis(self, a_nm, alpha_deg, beta_deg=0.0, convention="co_rotate"):
        b0 = self.rhombo_basis(a_nm, alpha_deg)
        if convention == "co_rotate":
            Ry, _ = self.rotation_parts(beta_deg)
            return b0 @ Ry.T
        if convention == "fixed_lattice":
            return b0
        raise ValueError("convention must be 'co_rotate' or 'fixed_lattice'")

    def site_indices(self, n=1):
        r = np.arange(-n, n + 1)
        return np.array(np.meshgrid(r, r, r, indexing="ij")).T.reshape(-1, 3)

    def supercell_sites(self, a_nm, alpha_deg, beta_deg=0.0, n=1,
                        convention="co_rotate", cutoff_nm=None):
        idx = self.site_indices(n)
        pos = idx @ self.lattice_basis(a_nm, alpha_deg, beta_deg, convention)
        dist_nm = np.linalg.norm(pos, axis=1) * 1e9
        keep = np.ones(len(idx), dtype=bool)
        if cutoff_nm is not None:
            keep = (dist_nm <= cutoff_nm) | (dist_nm < 1e-9)
        idx = idx[keep]
        pos = pos[keep]
        dist_nm = dist_nm[keep]

        center = np.where(np.all(idx == 0, axis=1))[0]
        if len(center) != 1:
            raise RuntimeError("Cluster must contain exactly one central site.")
        center = int(center[0])
        order = [center] + [i for i in np.argsort(dist_nm) if i != center]
        return idx[order], pos[order]

    def first_shell_count(self, a_nm, alpha_deg, n=1, tol_nm=0.03):
        idx, pos = self.supercell_sites(a_nm, alpha_deg, n=n, cutoff_nm=None)
        d_nm = np.linalg.norm(pos[1:], axis=1) * 1e9
        d0 = float(np.min(d_nm))
        return int(np.sum(np.abs(d_nm - d0) <= tol_nm)), d0

    def initial_orientations(self, idx, beta_deg=0.0, phi_deg=0.0):
        Ry, chain_hat = self.rotation_parts(beta_deg)
        Q0 = Ry @ self.Rbase
        Qs = np.empty((len(idx), 3, 3))
        layers = idx.sum(axis=1)
        for i, k in enumerate(layers):
            if abs(phi_deg) < 1e-14 or k == 0:
                Qs[i] = Q0
            else:
                Rphi = R.from_rotvec(np.radians(k * phi_deg) * chain_hat).as_matrix()
                Qs[i] = Rphi @ Q0
        return Qs

    def initial_state(self, a_nm, alpha_deg, beta_deg=0.0, phi_deg=0.0,
                      n=1, convention="co_rotate", cutoff_nm=None,
                      dipoles="field", dipole_rng=None):
        idx, pos = self.supercell_sites(a_nm, alpha_deg, beta_deg, n, convention, cutoff_nm)
        Qs = self.initial_orientations(idx, beta_deg, phi_deg)
        if dipoles == "field":
            mu = np.tile(self.field_hat, (len(idx), 1))
        elif dipoles == "easy":
            mu = np.array([Q @ self.e111 for Q in Qs])
        elif dipoles == "random":
            rr = rng if dipole_rng is None else dipole_rng
            mu = random_unit_vector(len(idx), rr)
        else:
            raise ValueError("dipoles must be field, easy, or random")
        return {"idx": idx, "pos": pos, "Q": Qs, "mu": mu}

    def support(self, Q, direction):
        d = np.asarray(direction, float)
        n = np.linalg.norm(d)
        if n < 1e-15:
            return 0.0
        dhat = d / n
        L1 = np.sum(np.abs(Q.T @ dhat))
        return (self.L * 0.5 - self.r_round) * L1 + self.r_round

    def surface_gap_pair(self, rij, Qi, Qj):
        rmag = np.linalg.norm(rij)
        if rmag < 1e-15:
            return -self.L
        rhat = rij / rmag
        return rmag - self.support(Qi, rhat) - self.support(Qj, -rhat)

    def anisotropy_energy(self, Qs, mu_hats, anis_model=None):
        model = self.anis_model if anis_model is None else anis_model
        if model == "none":
            return 0.0
        if model == "uniaxial111":
            easy = np.einsum("nij,j->ni", Qs, self.e111)
            dots = np.sum(mu_hats * easy, axis=1)
            return -self.K * self.vol * np.sum(dots**2)
        mu_body = np.array([Q.T @ mu for Q, mu in zip(Qs, mu_hats)])
        mx, my, mz = mu_body[:, 0], mu_body[:, 1], mu_body[:, 2]
        term = (mx * my) ** 2 + (mx * mz) ** 2 + (my * mz) ** 2
        if model == "cubic_first_raw":
            return -self.K * self.vol * np.sum(term)
        if model == "cubic_first_scaled":
            # Scaled so a <111> easy direction has the same depth as V6.02's
            # uniaxial -K V minimum. Raw cubic first order gives -K V / 3.
            return -3.0 * self.K * self.vol * np.sum(term)
        raise ValueError("anis_model must be none, uniaxial111, cubic_first_raw, or cubic_first_scaled")

    def energy_full_cluster(self, pos, Qs, mu_hats, a_nm=None, anis_model=None,
                            vdw_cut_factor=2.2, return_pairs=False):
        N = len(pos)
        Ez = -self.mu_mag * float(np.sum(mu_hats @ self.Bvec))
        Ea = self.anisotropy_energy(Qs, mu_hats, anis_model)
        Edd = 0.0
        Evdw = 0.0
        Ester = 0.0
        pair_rows = []
        a_m = None if a_nm is None else a_nm * 1e-9
        vox = None

        for i in range(N - 1):
            for j in range(i + 1, N):
                rij = pos[j] - pos[i]
                d = np.linalg.norm(rij)
                if d < 1e-15:
                    Ester += 1e25
                    continue
                rhat = rij / d
                mui = self.mu_mag * mu_hats[i]
                muj = self.mu_mag * mu_hats[j]
                Edd += self.Cmag * (
                    np.dot(mui, muj) / d**3
                    - 3.0 * np.dot(mui, rhat) * np.dot(muj, rhat) / d**3
                )

                gap = self.surface_gap_pair(rij, Qs[i], Qs[j])
                if gap < self.gap_m:
                    Ester += 0.5 * self.k_stiff * (self.gap_m - gap) ** 2

                use_vdw = a_m is None or d < vdw_cut_factor * a_m
                if use_vdw:
                    if vox is None:
                        vox = [self.vox_ref @ Q.T + p for Q, p in zip(Qs, pos)]
                    d2 = cdist(vox[i], vox[j], "sqeuclidean")
                    Evdw += self.scv * np.sum(1.0 / np.maximum(d2, 1e-19) ** 3)

                if return_pairs:
                    pair_rows.append((i, j, d * 1e9, gap * 1e9))

        out = {
            "Total": Ez + Ea + Edd + Evdw + Ester,
            "Zeeman": Ez,
            "Anisotropy": Ea,
            "Dipole": Edd,
            "VdW": Evdw,
            "Steric": Ester,
        }
        if return_pairs:
            out["pairs"] = pair_rows
        return out

    def single_energy(self, Q, mu_hat, anis_model=None):
        Ez = -self.mu_mag * float(np.dot(mu_hat, self.Bvec))
        Ea = self.anisotropy_energy(Q[None, :, :], mu_hat[None, :], anis_model)
        return Ez + Ea

    def pair_energy(self, pi, Qi, mui, pj, Qj, muj, a_nm=None, include_vdw=True):
        rij = pj - pi
        d = np.linalg.norm(rij)
        if d < 1e-15:
            return 1e25
        rhat = rij / d
        mui_v = self.mu_mag * mui
        muj_v = self.mu_mag * muj
        Edd = self.Cmag * (
            np.dot(mui_v, muj_v) / d**3
            - 3.0 * np.dot(mui_v, rhat) * np.dot(muj_v, rhat) / d**3
        )
        gap = self.surface_gap_pair(rij, Qi, Qj)
        Ester = 0.5 * self.k_stiff * (self.gap_m - gap) ** 2 if gap < self.gap_m else 0.0
        Evdw = 0.0
        if include_vdw:
            use_vdw = True if a_nm is None else d < 2.2 * a_nm * 1e-9
            if use_vdw:
                vi = self.vox_ref @ Qi.T + pi
                vj = self.vox_ref @ Qj.T + pj
                d2 = cdist(vi, vj, "sqeuclidean")
                Evdw = self.scv * np.sum(1.0 / np.maximum(d2, 1e-19) ** 3)
        return Edd + Ester + Evdw

    def local_energy(self, i, pos, Qs, mu_hats, a_nm=None, anis_model=None,
                     include_vdw=True):
        E = self.single_energy(Qs[i], mu_hats[i], anis_model)
        for j in range(len(pos)):
            if j == i:
                continue
            E += self.pair_energy(pos[i], Qs[i], mu_hats[i],
                                  pos[j], Qs[j], mu_hats[j],
                                  a_nm=a_nm, include_vdw=include_vdw)
        return E

    def magnetic_energy_only(self, Qs, mu_hats, anis_model=None):
        Ez = -self.mu_mag * float(np.sum(mu_hats @ self.Bvec))
        Ea = self.anisotropy_energy(Qs, mu_hats, anis_model)
        Edd = 0.0
        return Ez, Ea, Edd

    def magnetic_energy_for_state(self, pos, Qs, mu_hats, anis_model=None):
        Ez = -self.mu_mag * float(np.sum(mu_hats @ self.Bvec))
        Ea = self.anisotropy_energy(Qs, mu_hats, anis_model)
        Edd = 0.0
        N = len(pos)
        for i in range(N - 1):
            for j in range(i + 1, N):
                rij = pos[j] - pos[i]
                d = np.linalg.norm(rij)
                rhat = rij / d
                mui = self.mu_mag * mu_hats[i]
                muj = self.mu_mag * mu_hats[j]
                Edd += self.Cmag * (
                    np.dot(mui, muj) / d**3
                    - 3.0 * np.dot(mui, rhat) * np.dot(muj, rhat) / d**3
                )
        return {"Total": Ez + Ea + Edd, "Zeeman": Ez, "Anisotropy": Ea, "Dipole": Edd}

    def central_star_energy(self, a_nm, alpha_deg, beta_deg=0.0, phi_deg=0.0,
                            convention="co_rotate", moment_model="field_locked",
                            verbose=False):
        idx, pos = self.supercell_sites(a_nm, alpha_deg, beta_deg, n=1,
                                        convention=convention, cutoff_nm=None)
        Qs = self.initial_orientations(idx, beta_deg, phi_deg)
        center = 0
        c0 = pos[center]
        Q0 = Qs[center]
        easy0 = Q0 @ self.e111
        if moment_model == "field_locked":
            mu0_hat = self.field_hat
            mu_hats = np.tile(self.field_hat, (len(idx), 1))
        elif moment_model == "easy_axis_locked":
            mu_hats = np.array([Q @ self.e111 for Q in Qs])
            mu0_hat = mu_hats[center]
        else:
            raise ValueError("moment_model must be field_locked or easy_axis_locked")

        Ez = -self.mu_mag * float(np.dot(mu0_hat, self.Bvec))
        Ea = self.anisotropy_energy(Q0[None, :, :], mu0_hat[None, :], self.anis_model)
        Edd = 0.0
        Evdw = 0.0
        Ester = 0.0
        vox0 = self.vox_ref @ Q0.T + c0
        a_m = a_nm * 1e-9

        for j in range(1, len(idx)):
            rij = pos[j] - c0
            d = np.linalg.norm(rij)
            rhat = rij / d
            muj_hat = mu_hats[j]
            mui = self.mu_mag * mu0_hat
            muj = self.mu_mag * muj_hat
            Edd += 0.5 * self.Cmag * (
                np.dot(mui, muj) / d**3
                - 3.0 * np.dot(mui, rhat) * np.dot(muj, rhat) / d**3
            )
            gap = self.surface_gap_pair(rij, Q0, Qs[j])
            if gap < self.gap_m:
                Ester += 0.5 * self.k_stiff * (self.gap_m - gap) ** 2
            if d < 2.2 * a_m:
                voxj = self.vox_ref @ Qs[j].T + pos[j]
                d2 = cdist(vox0, voxj, "sqeuclidean")
                Evdw += 0.5 * self.scv * np.sum(1.0 / np.maximum(d2, 1e-19) ** 3)

        out = {
            "Total": Ez + Ea + Edd + Evdw + Ester,
            "Zeeman": Ez,
            "Anisotropy": Ea,
            "Dipole": Edd,
            "VdW": Evdw,
            "Steric": Ester,
        }
        return out if verbose else out["Total"]

    def optimise_v602_like(self, x0=None):
        if x0 is None:
            x0 = [self.L_nm * 2.0, 60.0]
        def obj(x):
            a_nm, alpha = float(x[0]), float(x[1])
            if a_nm < self.L_nm * 0.88 or not (45.0 < alpha < 135.0):
                return 1e25
            return self.central_star_energy(a_nm, alpha, 0.0, 0.0,
                                            convention="co_rotate",
                                            moment_model="field_locked")
        res = minimize(obj, x0, method="Nelder-Mead", tol=1e-6)
        return res.x, float(res.fun)

    def beta_angles(self, Qs):
        easy = np.einsum("nij,j->ni", Qs, self.e111)
        dots = np.clip(easy @ self.field_hat, -1.0, 1.0)
        return np.degrees(np.arccos(dots))

    def axial_tilt_deg(self, axis):
        axis = np.asarray(axis, float)
        norm = np.linalg.norm(axis)
        if norm < 1e-15:
            return np.nan
        dot = abs(float(np.dot(axis / norm, self.field_hat)))
        return float(np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))))

    def principal_axis_from_moment(self, moment):
        vals, vecs = np.linalg.eigh(moment)
        order = np.argsort(vals)[::-1]
        vals = vals[order]
        axis = vecs[:, order[0]]
        if np.dot(axis, self.field_hat) < 0.0:
            axis = -axis
        denom = float(np.sum(np.maximum(vals, 0.0)))
        orient_order = 0.0 if denom < 1e-30 else float((vals[0] - vals[1]) / denom)
        return axis, max(orient_order, 0.0), vals

    def sl_axis_pca(self, pos):
        centered = pos - np.mean(pos, axis=0, keepdims=True)
        moment = centered.T @ centered
        return self.principal_axis_from_moment(moment)

    def sl_tilt_metrics(self, pos):
        pca_axis, pca_order, pca_vals = self.sl_axis_pca(pos)
        return {
            "pca_axis": pca_axis,
            "pca_tilt_deg": self.axial_tilt_deg(pca_axis),
            "pca_order": pca_order,
            "pca_evals": pca_vals,
        }

    def body_axis_metrics(self, Qs):
        easy = np.einsum("nij,j->ni", Qs, self.e111)
        beta = self.beta_angles(Qs)
        mean_axis = np.mean(easy, axis=0)
        vector_order = float(np.linalg.norm(mean_axis))
        coherent_tilt = self.axial_tilt_deg(mean_axis)

        moment = easy.T @ easy
        nem_axis, nem_order, nem_vals = self.principal_axis_from_moment(moment)
        return {
            "local_beta_mean_deg": float(np.mean(beta)),
            "local_beta_std_deg": float(np.std(beta)),
            "coherent_axis": mean_axis / vector_order if vector_order > 1e-15 else np.full(3, np.nan),
            "coherent_tilt_deg": coherent_tilt,
            "vector_order": vector_order,
            "nematic_axis": nem_axis,
            "nematic_tilt_deg": self.axial_tilt_deg(nem_axis),
            "nematic_order": nem_order,
            "nematic_evals": nem_vals,
        }

    def axis_angle_deg(self, a, b, axial=True):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        na = np.linalg.norm(a)
        nb = np.linalg.norm(b)
        if na < 1e-15 or nb < 1e-15:
            return np.nan
        dot = float(np.dot(a / na, b / nb))
        if axial:
            dot = abs(dot)
        return float(np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))))

    def body_sl_mismatch_deg(self, pos, Qs):
        body = self.body_axis_metrics(Qs)
        sl = self.sl_tilt_metrics(pos)
        easy = np.einsum("nij,j->ni", Qs, self.e111)
        local_to_sl = [
            self.axis_angle_deg(e, sl["pca_axis"], axial=True)
            for e in easy
        ]
        return {
            "coherent_pca_deg": self.axis_angle_deg(
                body["coherent_axis"], sl["pca_axis"], axial=True
            ),
            "local_pca_mean_deg": float(np.mean(local_to_sl)),
            "local_pca_std_deg": float(np.std(local_to_sl)),
        }

    def surface_gap_stats(self, pos, Qs):
        gaps = []
        for i in range(len(pos) - 1):
            for j in range(i + 1, len(pos)):
                gaps.append(self.surface_gap_pair(pos[j] - pos[i], Qs[i], Qs[j]) * 1e9)
        gaps = np.asarray(gaps)
        return {
            "min_gap_nm": float(np.min(gaps)),
            "p05_gap_nm": float(np.percentile(gaps, 5.0)),
            "mean_gap_nm": float(np.mean(gaps)),
        }

    def layer_phi(self, idx, Qs):
        layers = idx.sum(axis=1)
        unique = np.array(sorted(set(layers)))
        angles = []
        used_layers = []
        for k in unique:
            mask = layers == k
            u_lab = np.einsum("nij,j->ni", Qs[mask], self.body_u)
            # azimuth around lab x; good diagnostic for co-rotated chains.
            ang = np.arctan2(u_lab[:, 2], u_lab[:, 1])
            mean_ang = np.angle(np.mean(np.exp(1j * ang)))
            angles.append(mean_ang)
            used_layers.append(k)
        if len(angles) < 2:
            return 0.0
        angles = np.unwrap(np.array(angles))
        used_layers = np.array(used_layers)
        dphi = np.diff(angles) / np.diff(used_layers)
        return float(np.degrees(np.mean(dphi)))


print("CollectiveNanocubeMC defined.")


def metropolis_accept(dE, kT, rng):
    return dE <= 0.0 or rng.random() < np.exp(-dE / kT)


def rotate_vector_small(v, max_angle_rad, rng):
    axis = random_unit_vector(1, rng)
    angle = rng.uniform(-max_angle_rad, max_angle_rad)
    out = R.from_rotvec(angle * axis).as_matrix() @ v
    return out / np.linalg.norm(out)


def rotate_body_small(Q, max_angle_deg, rng):
    axis = random_unit_vector(1, rng)
    angle = np.radians(rng.uniform(-max_angle_deg, max_angle_deg))
    return R.from_rotvec(angle * axis).as_matrix() @ Q


def rotate_cluster_small(pos, Qs, anchor, max_angle_deg, rng):
    axis = random_unit_vector(1, rng)
    angle = np.radians(rng.uniform(-max_angle_deg, max_angle_deg))
    Rg = R.from_rotvec(angle * axis).as_matrix()
    rel = pos - anchor
    new_pos = (Rg @ rel.T).T + anchor
    new_Qs = np.einsum("ab,nbc->nac", Rg, Qs)
    return new_pos, new_Qs


def rotate_bodies_about_coherent_axis(Qs, model, max_angle_deg, rng):
    body = model.body_axis_metrics(Qs)
    axis = body["coherent_axis"]
    if not np.all(np.isfinite(axis)):
        axis = model.field_hat
    angle = np.radians(rng.uniform(-max_angle_deg, max_angle_deg))
    Rg = R.from_rotvec(angle * axis).as_matrix()
    new_Qs = np.einsum("ab,nbc->nac", Rg, Qs)
    return new_Qs


def progress_range(n, enabled=False, desc="MC", backend="text", every=50):
    if not enabled:
        return range(n)
    if backend == "none":
        return range(n)
    if backend == "text":
        every = max(int(every), 1)

        def iterator():
            t0 = time.perf_counter()
            for c in range(n):
                yield c
                done = c + 1
                if done == 1 or done == n or done % every == 0:
                    elapsed = time.perf_counter() - t0
                    rate = done / max(elapsed, 1e-12)
                    eta = (n - done) / max(rate, 1e-12)
                    pct = 100.0 * done / max(n, 1)
                    print(f"{desc}: {done}/{n} ({pct:5.1f}%), "
                          f"elapsed {elapsed/60:.1f} min, ETA {eta/60:.1f} min",
                          flush=True)

        return iterator()
    try:
        from tqdm.auto import trange
        return trange(n, desc=desc, leave=True)
    except ImportError:
        print("tqdm is not installed; running without a progress bar.")
        return range(n)


def report_value(value):
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.6g}"
    return str(value)


def rows_to_text_table(rows, columns):
    if not rows:
        return ["(no rows)"]
    text_rows = [[report_value(row.get(col, "")) for col in columns] for row in rows]
    widths = [
        max(len(str(col)), *(len(row[i]) for row in text_rows))
        for i, col in enumerate(columns)
    ]
    header = " | ".join(str(col).ljust(widths[i]) for i, col in enumerate(columns))
    rule = "-+-".join("-" * w for w in widths)
    out = [header, rule]
    out.extend(" | ".join(row[i].ljust(widths[i]) for i in range(len(columns)))
               for row in text_rows)
    return out


def save_text_page(pdf, title, lines, fontsize=9, figsize=(8.5, 11.0)):
    fig = plt.figure(figsize=figsize)
    fig.text(0.05, 0.965, title, fontsize=14, weight="bold", va="top")
    fig.text(0.05, 0.925, "\n".join(lines), fontsize=fontsize,
             family="monospace", va="top")
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)


def base_report_lines():
    return [
        f"RUN_LABEL              = {RUN_LABEL}",
        f"RUN_PRESET             = {RUN_PRESET}",
        f"L_NM, A_NM, ALPHA_DEG  = {L_NM}, {A_NM}, {ALPHA_DEG}",
        f"BETA_INIT_DEG          = {BETA_INIT_DEG}",
        f"ANIS_MODEL             = {ANIS_MODEL}",
        f"DIPOLE_INIT_MODE       = {DIPOLE_INIT_MODE}",
        f"MOVE_POSITIONS         = {MOVE_POSITIONS}",
        f"MOVE_ORIENTATIONS      = {MOVE_ORIENTATIONS}",
        f"GLOBAL_COTILT          = {GLOBAL_COTILT_MOVES_PER_CYCLE} moves/cycle, "
        f"{GLOBAL_COTILT_STEP_DEG} deg",
        f"GLOBAL_GAMMA           = {GLOBAL_GAMMA_MOVES_PER_CYCLE} moves/cycle, "
        f"{GLOBAL_GAMMA_STEP_DEG} deg",
        f"N_MAG_PER_CYCLE        = {N_MAG_PER_CYCLE}",
        f"TRANS_STEP_NM          = {TRANS_STEP_NM}",
        f"ROT_STEP_DEG           = {ROT_STEP_DEG}",
        f"DIP_STEP_RAD           = {DIP_STEP_RAD}",
        f"B_FIELD_T              = {PHYS_500G['B_field']}",
        f"N_VOXEL                = {PHYS_500G['N_voxel']}",
    ]


def run_local_collective_mc(model, state, a_nm, n_cycles=500, n_equil=150,
                            n_mag_per_cycle=2, trans_step_nm=0.03,
                            rot_step_deg=2.0, dip_step_rad=0.30,
                            anis_model="cubic_first_raw", include_vdw=True,
                            move_positions=True, move_orientations=True,
                            n_cotilt_per_cycle=0, cotilt_step_deg=0.75,
                            n_gamma_per_cycle=0, gamma_step_deg=2.0,
                            show_progress=False, progress_label="MC",
                            progress_backend="text", progress_every=50,
                            rng=rng):
    """Local-delta MC for N=27 scale clusters."""
    idx = state["idx"].copy()
    pos = state["pos"].copy()
    Qs = state["Q"].copy()
    mu = state["mu"].copy()
    N = len(pos)
    center = 0
    E = model.energy_full_cluster(pos, Qs, mu, a_nm=a_nm, anis_model=anis_model)["Total"]

    traj_E = np.empty(n_cycles)
    traj_beta = np.empty(n_cycles)
    traj_body_tilt = np.empty(n_cycles)
    traj_body_order = np.empty(n_cycles)
    traj_body_sl_mismatch = np.empty(n_cycles)
    traj_sl_pca_tilt = np.empty(n_cycles)
    traj_sl_pca_order = np.empty(n_cycles)
    traj_gap_min = np.empty(n_cycles)
    traj_gap_p05 = np.empty(n_cycles)
    traj_muB = np.empty(n_cycles)
    acc_m = try_m = acc_c = try_c = acc_g = try_g = acc_d = try_d = 0

    for c in progress_range(n_cycles, enabled=show_progress, desc=progress_label,
                            backend=progress_backend, every=progress_every):
        for i in rng.permutation(np.arange(N)):
            if i == center:
                continue
            old_pos = pos[i].copy()
            old_Q = Qs[i].copy()
            old_local = model.local_energy(i, pos, Qs, mu, a_nm=a_nm,
                                           anis_model=anis_model, include_vdw=include_vdw)
            if move_positions:
                pos[i] = pos[i] + random_unit_vector(1, rng) * rng.uniform(0, trans_step_nm) * 1e-9
            if move_orientations:
                Qs[i] = rotate_body_small(Qs[i], rot_step_deg, rng)
            new_local = model.local_energy(i, pos, Qs, mu, a_nm=a_nm,
                                           anis_model=anis_model, include_vdw=include_vdw)
            dE = new_local - old_local
            if metropolis_accept(dE, model.kT, rng):
                E += dE
                acc_m += 1
            else:
                pos[i] = old_pos
                Qs[i] = old_Q
            try_m += 1

        for _ in range(n_cotilt_per_cycle):
            old_pos = pos.copy()
            old_Qs = Qs.copy()
            new_pos, new_Qs = rotate_cluster_small(
                pos, Qs, pos[center].copy(), cotilt_step_deg, rng
            )
            new_E = model.energy_full_cluster(new_pos, new_Qs, mu, a_nm=a_nm,
                                              anis_model=anis_model)["Total"]
            dE = new_E - E
            if metropolis_accept(dE, model.kT, rng):
                pos = new_pos
                Qs = new_Qs
                E = new_E
                acc_c += 1
            else:
                pos = old_pos
                Qs = old_Qs
            try_c += 1

        for _ in range(n_gamma_per_cycle):
            old_Qs = Qs.copy()
            new_Qs = rotate_bodies_about_coherent_axis(Qs, model, gamma_step_deg, rng)
            new_E = model.energy_full_cluster(pos, new_Qs, mu, a_nm=a_nm,
                                              anis_model=anis_model)["Total"]
            dE = new_E - E
            if metropolis_accept(dE, model.kT, rng):
                Qs = new_Qs
                E = new_E
                acc_g += 1
            else:
                Qs = old_Qs
            try_g += 1

        for _ in range(n_mag_per_cycle * N):
            i = int(rng.integers(N))
            old = mu[i].copy()
            old_local = model.local_energy(i, pos, Qs, mu, a_nm=a_nm,
                                           anis_model=anis_model, include_vdw=False)
            mu[i] = rotate_vector_small(mu[i], dip_step_rad, rng)
            new_local = model.local_energy(i, pos, Qs, mu, a_nm=a_nm,
                                           anis_model=anis_model, include_vdw=False)
            dE = new_local - old_local
            if metropolis_accept(dE, model.kT, rng):
                E += dE
                acc_d += 1
            else:
                mu[i] = old
            try_d += 1

        traj_E[c] = E / model.kT
        traj_beta[c] = np.mean(model.beta_angles(Qs))
        body = model.body_axis_metrics(Qs)
        sl = model.sl_tilt_metrics(pos)
        mismatch = model.body_sl_mismatch_deg(pos, Qs)
        gaps = model.surface_gap_stats(pos, Qs)
        traj_body_tilt[c] = body["coherent_tilt_deg"]
        traj_body_order[c] = body["vector_order"]
        traj_body_sl_mismatch[c] = mismatch["coherent_pca_deg"]
        traj_sl_pca_tilt[c] = sl["pca_tilt_deg"]
        traj_sl_pca_order[c] = sl["pca_order"]
        traj_gap_min[c] = gaps["min_gap_nm"]
        traj_gap_p05[c] = gaps["p05_gap_nm"]
        traj_muB[c] = np.mean(np.abs(mu @ model.Bhat))

    eq = slice(n_equil, None)
    return {
        "idx": idx, "pos": pos, "Q": Qs, "mu": mu,
        "E_mean": float(np.mean(traj_E[eq])),
        "beta_mean": float(np.mean(traj_beta[eq])),
        "body_tilt_mean": float(np.mean(traj_body_tilt[eq])),
        "body_order_mean": float(np.mean(traj_body_order[eq])),
        "body_sl_mismatch_mean": float(np.mean(traj_body_sl_mismatch[eq])),
        "sl_pca_tilt_mean": float(np.mean(traj_sl_pca_tilt[eq])),
        "sl_pca_order_mean": float(np.mean(traj_sl_pca_order[eq])),
        "gap_min_mean": float(np.mean(traj_gap_min[eq])),
        "gap_p05_mean": float(np.mean(traj_gap_p05[eq])),
        "muB_mean": float(np.mean(traj_muB[eq])),
        "acc_mech": acc_m / max(try_m, 1),
        "acc_cotilt": acc_c / max(try_c, 1),
        "acc_gamma": acc_g / max(try_g, 1),
        "acc_dip": acc_d / max(try_d, 1),
        "traj_E": traj_E, "traj_beta": traj_beta,
        "traj_body_tilt": traj_body_tilt,
        "traj_body_order": traj_body_order,
        "traj_body_sl_mismatch": traj_body_sl_mismatch,
        "traj_sl_pca_tilt": traj_sl_pca_tilt,
        "traj_sl_pca_order": traj_sl_pca_order,
        "traj_gap_min_nm": traj_gap_min,
        "traj_gap_p05_nm": traj_gap_p05,
        "traj_muB": traj_muB,
    }


print("V9.03 local-delta MC helper defined.")
