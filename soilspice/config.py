"""Paths and global settings.

This is the only file with machine-specific paths. Everything else reads from here.

Inputs
    data/netlists/       LTspice netlists of the five modules (.net)
    data/climate/        ERA5 hourly series (Open-Meteo) and MeteoSwiss station NEU
    data/field_test/     IV-tracer export of the roof field test
    DATA_DIR             flash-test campaign (.asc files), kept outside the repository (300 MB)

Outputs, one folder per chapter of the paper (paper figures at the top of each folder, the data behind them
and the diagnostic plots in data/):
    results/1_validation/      chapter 1, experimental validation of the model
    results/2_simulation/      chapter 2, electrical resilience across modules
    results/3_field_test/      chapter 3, roof field test
    results/4_economics/       chapter 4, energy yield and cleaning costs
"""
from pathlib import Path

# --- machine-specific ------------------------------------------------------------------------------------------
DATA_DIR = Path("/Users/lorrainenaux/Desktop/Valorisation-master/Data/IV-curve")   # flash-test .asc files
LTSPICE_EXE = Path("/Applications/LTspice.app/Contents/MacOS/LTspice")

# --- inside the repository --------------------------------------------------------------------------------------
REPO_DIR = Path(__file__).resolve().parent.parent
INPUT_DIR = REPO_DIR / "data"
NETLIST_DIR = INPUT_DIR / "netlists"
CLIMATE_DIR = INPUT_DIR / "climate"
FIELD_CSV = INPUT_DIR / "field_test" / "rooftop_SPR_2026-09-21_2026-10-15.csv"

OUTPUT_DIR = REPO_DIR / "results"
CH1_DIR = OUTPUT_DIR / "1_validation"
CH2_DIR = OUTPUT_DIR / "2_simulation"
CH3_DIR = OUTPUT_DIR / "3_field_test"
CH4_DIR = OUTPUT_DIR / "4_economics"

# chapter 1
VALIDATION_DIR = CH1_DIR / "data" / "validation"          # per module: <MODULE>/
CROSS_MODULE_DIR = VALIDATION_DIR / "cross_module"
UNCERTAINTY_DIR = CH1_DIR / "data" / "uncertainty"
# chapter 2
SWEEP_DIR = CH2_DIR / "data"                               # loss maps at equal A_cov, resilience metrics
EQUAL_HEIGHT_DIR = SWEEP_DIR / "equal_height"              # losses at equal band height (used by chapters 3, 4)
LAYOUT_DIR = SWEEP_DIR / "layout"                          # drawings of the bands on each module
EQUAL_HEIGHT_CSV = EQUAL_HEIGHT_DIR / "equal_height_sweep.csv"
# chapter 3
FIELD_DIR = CH3_DIR / "data"
# chapter 4
ECON_DIR = CH4_DIR / "data"

# Patched netlists and .raw files. Can be deleted at any time.
RUNS_DIR = REPO_DIR / "spice_runs"

# --- LTspice conventions shared by all module netlists ----------------------------------------------------------
VOLTAGE_TRACE = "v1"       # swept source voltage = module voltage
CURRENT_TRACE = "I(V1)"
CURRENT_SIGN = 1           # set to -1 if the current comes out negative
TEMPERATURE_PARAM = "T"    # name of the `.param T=...` line in each netlist
LTSPICE_TIMEOUT_S = 600

# --- measurement handling ---------------------------------------------------------------------------------------
# Shaded flashes are corrected to the temperature of the clean flash (read in the .asc file) before being compared
# with the clean flash and with SPICE.
NORMALIZE_TO_REF_TEMPERATURE = True
