"""Description of the PV modules: geometry, LTspice netlist and cell numbering.

A module is described once, in portrait orientation, with rows counted from
the bottom edge. Landscape is derived from it in `geometry.py`.

To add a module: write its cell-map function, add a `Module(...)` entry to
MODULES, and check the result with `scripts/check_layout.py` before running
any simulation.
"""
from dataclasses import dataclass
from functools import cached_property
from typing import Callable, Optional, Union


@dataclass(frozen=True)
class PVsystShuntLaw:
    """Cell shunt resistance at STC and the prescribed PVsyst exponential shape."""
    r_stc: float
    zero_ratio: float = 4.0
    exponent: float = 5.5

from . import config


@dataclass(frozen=True)
class Module:
    # identity (used in legends)
    name: str
    pretty_name: str
    technology: str

    # portrait geometry [mm]
    n_cols: int
    n_rows: int
    width: float
    height: float
    cell_width: float
    cell_height: float
    gap: float = 0.0
    center_gap: float = 0.0            # extra gap between the two half-modules (between rows in portrait)
    margins: tuple = (0.0, 0.0, 0.0, 0.0)  # glass edge -> cells: (left, right, bottom, top)
    frame: tuple = (0.0, 0.0, 0.0, 0.0)    # frame width: (left, right, bottom, top)

    # electrical model
    netlist: str = ""                  # file name inside config.NETLIST_DIR
    subckt: str = ""                   # name of the cell sub-circuit in the netlist
    cell_prefix: str = "XX"            # instance names look like XX12
    cell_map: Callable[[], dict] = None   # {(row, col): instance}, see row0_is_bottom;
                                          # instance None = inactive cell (drawn, not simulated)
    row0_is_bottom: bool = True        # False if the cell map counts rows from the top
    bypass_groups: tuple = ()          # (first, last) portrait column of each bypass sub-string
    bypass_axis: str = "cols"          # "rows" if the groups are ranges of portrait rows
    junction_box: Optional[str] = None # "top" if the j-box sits on the top edge (drawing only)

    # shading parameters written on each cell line:
    #   "fcov_tau"   -> f_cov=<covered fraction> tau_I=<transmittance>   (Bishop cells)
    #   "fcov_only"  -> f_cov=<covered fraction x (1 - tau_I)>          (cells with a single opacity factor)
    shading_params: str = "fcov_tau"

    # reference conditions
    p_clean: float = 0.0               # SPICE clean power [W], only used for the power ratio
    t_ref: float = 25.0                # SPICE temperature [degC] that reproduces the clean flash (fit value)
    temperature: str = "param"         # how t_ref is written: "param" -> .param T=..., "directive" -> .temp ...
    kv: float = -0.0031                # relative temperature coefficient of Voc [1/K]
    kc: float = 0.0004                 # relative temperature coefficient of Isc [1/K]

    # Optional PVsyst law on effective cell irradiance, or the legacy power-law
    # tuple (Rp_ref [Ohm], G_ref [W/m2], m). None keeps the template's fixed Rp.
    rp_law: Optional[Union[tuple, PVsystShuntLaw]] = None


    @property
    def label(self):
        return f"{self.pretty_name} ({self.technology})"

    @property
    def netlist_path(self):
        return config.NETLIST_DIR / self.netlist

    @property
    def output_dir(self):
        """Validation results of this module."""
        return config.VALIDATION_DIR / self.name

    @property
    def offsets(self):
        """Distance from the module edge to the first cell, margin + frame."""
        return tuple(m + f for m, f in zip(self.margins, self.frame))

    @cached_property
    def instances(self):
        """{(row, col): SPICE instance} of the active cells."""
        return {k: v for k, v in self.cell_map().items() if v is not None}

    @cached_property
    def inactive_cells(self):
        """Portrait (row, col) positions that hold a cell strip that produces no
        current (matrix shingle contact rows): covered by the band like the
        others, but not in the netlist."""
        return frozenset(k for k, v in self.cell_map().items() if v is None)

    def _key(self, row, col):
        return (self.n_rows - 1 - row, col) if not self.row0_is_bottom else (row, col)

    def instance(self, row, col):
        """SPICE instance of the portrait cell (row, col), row 0 = bottom row."""
        return self.instances.get(self._key(row, col))

    def is_inactive(self, row, col):
        return self._key(row, col) in self.inactive_cells


# --- cell maps -----------------------------------------------------------------
# Each function returns {(row, col): instance name} following the numbering of
# the LTspice schematic of that module. They were checked one by one against
# the schematics.

def yingli_cells():
    column_base = [11, 1, 31, 21, 51, 41]
    return {(row, col): f"XX{base + 9 - row}"
            for col, base in enumerate(column_base) for row in range(10)}


def jasolar_cells():
    # 6 columns x 18 rows, two half-modules of 9 rows; columns are wired in pairs.
    cells = {}
    for col in range(6):
        base = (col // 2) * 18
        left_of_pair = col % 2 == 0
        for k in range(9):
            row_top = 9 + k
            if left_of_pair:
                bottom, top = base + 1 + k, base + 72 - row_top
            else:
                bottom, top = base + 18 - k, base + 81 - row_top
            cells[(k, col)] = f"XX{bottom}"
            cells[(row_top, col)] = f"XX{top}"
    return cells


def longi_cells():
    # same numbering per bypass diode as JA Solar 
    return jasolar_cells()


def sunpower_cells():
    column_base = [1, 24, 25, 48, 49, 72, 73, 96]
    cells = {}
    for col, base in enumerate(column_base):
        step = 1 if col % 2 == 0 else -1      # serpentine wiring
        for row in range(12):
            cells[(row, col)] = f"XX{base + step * row}"
    return cells


def matrix_cells():
    # Installed (landscape) module: 2 blocks side by side, one bypass diode per
    # block. Each block has 14 shingled rows of 9 half-cells: 13 active rows and
    # one contact row on the + busbar (front side not connected, it only
    # conducts -> inactive, no current). The contact row is at the TOP of the
    # left block (A) and at the BOTTOM of the right block (B); turning the module
    # by 180 deg only swaps A and B.
    # Netlist names XXL<row>_<k> (left block, row 1 at the bottom) and
    # XXR<row>_<k> (right block, row 13 just above the contact row), k = 0..8
    # from left to right in each block.
    # Portrait = the landscape module turned so that its left edge is at the
    # bottom: portrait row i = landscape column i, portrait column j = landscape row j.
    cells = {}
    for i in range(18):                      # landscape column, left to right
        for j in range(14):                  # landscape row position, from the bottom
            if i < 9:                        # block A: rows 1..13, contact row on top
                cells[(i, j)] = f"XXL{j + 1}_{i}" if j < 13 else None
            else:                            # block B: contact row at the bottom, rows 13..1
                cells[(i, j)] = f"XXR{14 - j}_{i - 9}" if j > 0 else None
    return cells


# --- the modules ------------------------------------------------------------------

MODULES = {
    "LONGI": Module(
        name="LONGI", pretty_name="LONGi", technology="HJT BC",
        n_cols=6, n_rows=18, width=1134.0, height=1800.0,
        cell_width=182.2, cell_height=95.8, gap=1.0, center_gap=13.0,
        margins=(2.0, 2.0, 14.0, 12.0), frame=(10.0, 10.0, 10.0, 10.0),
        netlist="Longi-Bishop.net", subckt="pvcell-longi", cell_map=longi_cells,
        bypass_groups=((0, 1), (2, 3), (4, 5)),
        p_clean=479.082, t_ref=25.0, kv=-0.00210, kc=0.0005,
    ),
    "SUNPOWER": Module(
        name="SUNPOWER", pretty_name="SunPower", technology="IBC",
        n_cols=8, n_rows=12, width=1046.0, height=1558.0,
        cell_width=125.0, cell_height=125.0, gap=2.0,
        margins=(5.0, 5.0, 7.0, 7.0), frame=(10.0, 10.0, 10.0, 10.0),
        netlist="SUNPOWER-BISHOP-NEW.net", subckt="pvcell-sunpower-bishop",
        cell_map=sunpower_cells, row0_is_bottom=False,
        bypass_groups=((0, 1), (2, 5), (6, 7)), junction_box="top",
        p_clean=341.191, t_ref=25.4, kv=-0.002455, kc=0.000548,
    ),
    "YINGLI": Module(
        name="YINGLI", pretty_name="Yingli", technology="n-Si",
        n_cols=6, n_rows=10, width=990.0, height=1650.0,
        cell_width=156.0, cell_height=156.0, gap=2.0,
        margins=(6.0, 6.0, 17.0, 25.0), frame=(13.0, 13.0, 12.0, 13.0),
        netlist="Yingli_SOLAR_Bishop.net", subckt="pvcell_bishop", cell_map=yingli_cells,
        bypass_groups=((0, 1), (2, 3), (4, 5)), junction_box="top",
        p_clean=263.392, t_ref=24.8,
        # kv, kc: datasheet values not confirmed -> generic defaults
        rp_law=None,  # original constant cell shunt resistance: 5.99 ohm in the netlist
    ),
    "JASOLAR": Module(
        name="JASOLAR", pretty_name="JA Solar", technology="TOPCon",
        n_cols=6, n_rows=18, width=1134.0, height=1762.0,
        cell_width=181.85, cell_height=94.17, gap=1.0, center_gap=6.0,
        margins=(5.0, 5.0, 11.0, 13.0), frame=(10.0, 10.0, 10.0, 10.0),
        netlist="JASolar-bishop2.net", subckt="pvcell_jasolar", cell_map=jasolar_cells,
        bypass_groups=((0, 1), (2, 3), (4, 5)),
        p_clean=413.646, t_ref=24.2, kv=-0.0031, kc=0.0004,
    ),

}

# Plotting order and identity colours (the same in every figure).
# Fraunhofer ISE matrix shingle (noise-barrier module 250240ISE0031/32), Klasen et al. 2022.
# Installed with the long side horizontal (drawing 1980 x 386 mm) = landscape here.
# Cells: visible pitch of the shingled strips, (1980 - 2x33 - 20) / 18 and (386 - 2x16.95) / 14
# (13 active rows + 1 inactive contact row per block, assumed of the same height).
# No flash-test campaign: simulation only (validation skips it).
MODULES["MATRIX"] = Module(
    name="MATRIX", pretty_name="Matrix shingle", technology="Fraunhofer ISE",
    n_cols=14, n_rows=18, width=386.0, height=1980.0,
    cell_width=(386.0 - 2 * 16.95) / 14, cell_height=(1980.0 - 2 * 33.0 - 20.0) / 18,
    gap=0.0, center_gap=20.0,
    margins=(16.95, 16.95, 33.0, 33.0),      # (left, right, bottom, top) in portrait
    netlist="matrix-shingle-module.net", subckt="cell-matrix-shingle", cell_map=matrix_cells,
    bypass_groups=((0, 8), (9, 17)), bypass_axis="rows",
    shading_params="fcov_only", temperature="directive", t_ref=25.0,
)

MODULE_ORDER = ["LONGI", "SUNPOWER", "YINGLI", "JASOLAR", "MATRIX"]


def get(name):
    try:
        return MODULES[name.upper()]
    except KeyError:
        raise KeyError(f"Unknown module {name!r}. Known: {', '.join(MODULES)}") from None
