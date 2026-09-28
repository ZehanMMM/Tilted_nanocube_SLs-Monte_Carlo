"""Worst-case vdW attraction reachable at the steric wall with sharp cores >= D0 apart."""
import numpy as np, sys, time
from scipy.optimize import brentq
from scipy.spatial.transform import Rotation
sys.path.insert(0, '.')
from v904.vdw import HamakerEvaluator
from v904 import geometry
L = 16e-9; kT = 1.380649e-23 * 298.15; r_round = 1.5e-9; gap0 = 1.8e-9
ev = HamakerEvaluator()
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
D0 = float(sys.argv[2]) * 1e-9 if len(sys.argv) > 2 else 0.165e-9
best = []
t0 = time.time()
for trial in range(int(sys.argv[3]) if len(sys.argv) > 3 else 400):
    QA = np.eye(3)
    QB = Rotation.random(random_state=rng).as_matrix()
    u = rng.normal(size=3); u /= np.linalg.norm(u)
    # place B on the steric wall along u
    f = lambda s: geometry.support_gap(s * u, QA, QB, L, r_round)[0] - gap0
    s_wall = brentq(f, 0.5 * L, 3 * L)
    s = s_wall
    # push out until sharp cores are >= D0 apart
    for _ in range(60):
        if not geometry.cores_overlap(s * u, QA, QB, L, tol_m=0.0)[0]:
            d = geometry.core_distance(np.zeros(3), QA, s * u, QB, L)
            if d >= D0:
                break
        s += 0.02e-9
    e = ev.pair_energy(np.zeros(3), QA, s * u, QB) / kT
    best.append((e, s * 1e9, geometry.core_distance(np.zeros(3), QA, s * u, QB, L) * 1e9))
best.sort()
print(f'{len(best)} trials in {time.time() - t0:.0f}s; D0 = {D0 * 1e9} nm')
for b in best[:8]:
    print('E %.3f kBT  centre %.2f nm  core %.3f nm' % b)
print('fraction of wall contacts with core < 0.5 nm at the wall:',
      np.mean([b[2] < 0.5 for b in best]))
