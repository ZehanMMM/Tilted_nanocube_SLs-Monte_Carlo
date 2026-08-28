model = CollectiveNanocubeMC(L_NM, PHYS_500G, anis_model=ANIS_MODEL)


def make_initial_state(seed=INIT_SEED):
    return model.initial_state(A_NM, ALPHA_DEG, beta_deg=BETA_INIT_DEG,
                               phi_deg=PHI_INIT_DEG, n=1,
                               convention="co_rotate",
                               cutoff_nm=None, dipoles=DIPOLE_INIT_MODE,
                               dipole_rng=np.random.default_rng(seed))


state0 = make_initial_state(INIT_SEED)

E0 = model.energy_full_cluster(state0["pos"], state0["Q"], state0["mu"],
                               a_nm=A_NM, anis_model=ANIS_MODEL,
                               return_pairs=True)
gaps = np.array([p[3] for p in E0["pairs"]])
dist = np.linalg.norm(state0["pos"] - state0["pos"][0], axis=1) * 1e9
first = np.min(dist[1:])
first_count = int(np.sum(np.abs(dist[1:] - first) < 0.04))
beta0 = float(np.mean(model.beta_angles(state0["Q"])))
body0 = model.body_axis_metrics(state0["Q"])
sl0 = model.sl_tilt_metrics(state0["pos"])
mismatch0 = model.body_sl_mismatch_deg(state0["pos"], state0["Q"])
gap0 = model.surface_gap_stats(state0["pos"], state0["Q"])

print("Initial N27 cluster")
print(f"  initial seed = {INIT_SEED}, dipoles = {DIPOLE_INIT_MODE}")
print(f"  N = {len(state0['pos'])}")
print(f"  first shell around center = {first_count}, first distance = {first:.3f} nm")
print(f"  min pair gap = {np.min(gaps):.3f} nm")
print(f"  local cube-body <beta> = {beta0:.3f} deg")
print(f"  coherent cube-body tilt = {body0['coherent_tilt_deg']:.3f} deg, "
      f"body order = {body0['vector_order']:.3f}")
print(f"  SL tilt PCA = {sl0['pca_tilt_deg']:.3f} deg, order = {sl0['pca_order']:.3f}")
print(f"  coherent body-SL mismatch = {mismatch0['coherent_pca_deg']:.3f} deg")
print(f"  surface gap min/p05 = {gap0['min_gap_nm']:.3f} / {gap0['p05_gap_nm']:.3f} nm")
for key in ["Zeeman", "Anisotropy", "Dipole", "VdW", "Steric", "Total"]:
    print(f"  {key:<12} {E0[key]/model.kT:+12.3f} kBT")

assert len(state0["pos"]) == 27
assert np.min(gaps) > 0.0
