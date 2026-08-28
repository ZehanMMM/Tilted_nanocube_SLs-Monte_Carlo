import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import cdist
import scipy.constants as const

# =============================================================================
# User-editable parameters
# =============================================================================
DEFAULT_SEED = 1
INIT_SEED = DEFAULT_SEED
BENCH_SEED = DEFAULT_SEED + 1000
FULL_RUN_SEED = DEFAULT_SEED
N_CYCLES = 2000
N_EQUIL = 500
RUN_MULTI_SEED_SCAN = False
MULTI_SEEDS = [DEFAULT_SEED]  # e.g. [1, 2, 3]
MULTI_CYCLES = N_CYCLES
MULTI_EQUIL = N_EQUIL
MULTI_RESEED_INITIAL_STATE = True
MOVE_POSITIONS = True
MOVE_ORIENTATIONS = True
SHOW_PROGRESS_BARS = True
SHOW_BENCH_PROGRESS = False
PROGRESS_BACKEND = "text"  # text, tqdm, none
PROGRESS_EVERY_CYCLES = 50
SAVE_REPORT_PDFS = True
SHOW_FIGURES = True
FULL_RUN_FIG_PDF = "V903_Fig1_N27_2000cycles.pdf"
FULL_RUN_REPORT_PDF = "V903_FullRun_Report.pdf"
MULTI_SEED_REPORT_PDF = "V903_MultiSeed_Report.pdf"

EXP_DATA = {
    16.0: {"d": 21.0, "alpha": 74.2},
    18.0: {"d": 22.2, "alpha": 73.0},
    19.5: {"d": 29.0, "alpha": 60.0},
    29.6: {"d": 39.6, "alpha": 60.0},
}

DEFAULT_SIZE_NM = 16.0
RUN_PRESET = "size_default"  # size_default, 19p5_tilted_candidate, 29p6_tilted_candidate
DEFAULT_BETA_INIT_DEG = 0.0
DEFAULT_PHI_INIT_DEG = 0.0

TILTED_PRESETS = {
    "19p5_tilted_candidate": {
        "L_nm": 19.5,
        "a_nm": 31.0,
        "alpha_deg": 60.0,
        "beta_init_deg": 20.0,
        "phi_init_deg": 0.0,
        "label": "L=19.5, a=31.0, alpha=60, beta_init=20",
    },
    "29p6_tilted_candidate": {
        "L_nm": 29.6,
        "a_nm": EXP_DATA[29.6]["d"],
        "alpha_deg": 60.0,
        "beta_init_deg": 20.0,
        "phi_init_deg": 0.0,
        "label": "L=29.6, a=39.6, alpha=60, beta_init=20",
    },
}


def preset_from_size(size_nm, beta_init_deg=0.0, phi_init_deg=0.0):
    size_nm = float(size_nm)
    if size_nm not in EXP_DATA:
        raise KeyError(f"DEFAULT_SIZE_NM={size_nm} is not in EXP_DATA")
    exp = EXP_DATA[size_nm]
    return {
        "L_nm": size_nm,
        "a_nm": exp["d"],
        "alpha_deg": exp["alpha"],
        "beta_init_deg": float(beta_init_deg),
        "phi_init_deg": float(phi_init_deg),
        "label": f"L={size_nm:g}, a={exp['d']}, alpha={exp['alpha']}, beta_init={beta_init_deg}",
    }


def resolve_run_config(run_preset=None, default_size_nm=None):
    if run_preset is None:
        run_preset = RUN_PRESET
    if default_size_nm is None:
        default_size_nm = DEFAULT_SIZE_NM
    if run_preset == "size_default":
        return preset_from_size(default_size_nm, DEFAULT_BETA_INIT_DEG, DEFAULT_PHI_INIT_DEG)
    if run_preset in TILTED_PRESETS:
        return dict(TILTED_PRESETS[run_preset])
    allowed = ["size_default", *TILTED_PRESETS]
    raise KeyError(f"Unknown RUN_PRESET={run_preset!r}; choose one of {allowed}")

PHYS_V602 = {
    "Ms": 285000,
    "B_field": 0.5,      # T = 5000 G, original V6.02
    "K_ani": 2.0e4,
    "Hamaker": 2.0e-20,
    "N_voxel": 4,
    "gap_nm": 1.8,
    "k_stiff": 1.0e8,
    "roundness_nm": 1.5,
}

PHYS_500G = {**PHYS_V602, "B_field": 0.05}  # 500 G = 0.05 T

T_K = 298.15
RUN_CONFIG = resolve_run_config()
L_NM = RUN_CONFIG["L_nm"]
A_NM = RUN_CONFIG["a_nm"]
ALPHA_DEG = RUN_CONFIG["alpha_deg"]
BETA_INIT_DEG = RUN_CONFIG["beta_init_deg"]
PHI_INIT_DEG = RUN_CONFIG["phi_init_deg"]
RUN_LABEL = RUN_CONFIG["label"]
ANIS_MODEL = "cubic_first_raw"
DIPOLE_INIT_MODE = "random"  # random, field, easy

N_MAG_PER_CYCLE = 1
TRANS_STEP_NM = 0.03
ROT_STEP_DEG = 2.0
DIP_STEP_RAD = 0.30
GLOBAL_COTILT_MOVES_PER_CYCLE = 1
GLOBAL_COTILT_STEP_DEG = 0.75
GLOBAL_GAMMA_MOVES_PER_CYCLE = 1
GLOBAL_GAMMA_STEP_DEG = 2.0

rng = np.random.default_rng(INIT_SEED)

print("V9.03 N27 collective MC loaded.")
print(f"B = {PHYS_500G['B_field']} T = 500 G")
print(f"default seed={DEFAULT_SEED}; init={INIT_SEED}, bench={BENCH_SEED}, full={FULL_RUN_SEED}")
print(f"preset={RUN_PRESET}: {RUN_LABEL}")
print(f"L={L_NM} nm, a={A_NM} nm, alpha={ALPHA_DEG} deg, beta_init={BETA_INIT_DEG} deg")
print(f"cycles={N_CYCLES}, equil={N_EQUIL}, anis_model={ANIS_MODEL}")
print(f"global co-tilt moves/cycle={GLOBAL_COTILT_MOVES_PER_CYCLE}, step={GLOBAL_COTILT_STEP_DEG} deg")
print(f"global gamma moves/cycle={GLOBAL_GAMMA_MOVES_PER_CYCLE}, step={GLOBAL_GAMMA_STEP_DEG} deg")
