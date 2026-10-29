"""Running LTspice on a module netlist with a given shading pattern.

For each simulation the template netlist (exported from the LTspice schematic,
View > SPICE Netlist) is copied, every cell line gets its own f_cov / tau_I,
the temperature is set, and LTspice is run in batch mode. The I-V curve is
read back from the .raw file.
"""
import re
import subprocess
from dataclasses import dataclass

import numpy as np

from . import config
from .modules import PVsystShuntLaw

_PARAM_TOKEN = re.compile(r"(\w+)\s*=\s*(\S+)")
IV_KEYS = ("Isc_A", "Voc_V", "Imp_A", "Vmp_V", "Pmax_W", "FF_percent")


# --- netlist text ----------------------------------------------------------------

def read_netlist(path):
    """LTspice writes netlists in UTF-16 on some systems, UTF-8 on others."""
    raw = path.read_bytes()
    utf16 = raw[:2] in (b"\xff\xfe", b"\xfe\xff") or b"\x00" in raw[:200]
    return raw.decode("utf-16" if utf16 else "utf-8", errors="replace")


def _cell_line_regex(module):
    return re.compile(
        rf"^(?P<ref>{re.escape(module.cell_prefix)}\S*)\s+(?P<n1>\S+)\s+(?P<n2>\S+)\s+"
        rf"(?P<subckt>{re.escape(module.subckt)})\b(?P<rest>.*)$", re.IGNORECASE)


def shaded_rp(module, f_cov, tau_i, min_tau=1e-3):
    """Cell shunt resistance; PVsyst uses the same area-weighted light as Iph.

    The legacy optional power law is retained for existing callers.
    """
    if isinstance(module.rp_law, PVsystShuntLaw):
        law = module.rp_law
        g = 1.0 - f_cov * (1.0 - tau_i)
        if not 0.0 <= g <= 1.0:
            raise ValueError(f"relative cell irradiance outside [0, 1]: {g}")
        r_zero = law.zero_ratio * law.r_stc
        decay = np.exp(-law.exponent)
        base = (law.r_stc - r_zero * decay) / (1.0 - decay)
        return float(base + (r_zero - base) * np.exp(-law.exponent * g))
    rp_ref, _g_ref, m = module.rp_law
    if f_cov <= 0:
        return rp_ref
    return rp_ref / max(tau_i, min_tau) ** m


def patch_cells(text, module, cell_fcov, tau_i):
    """Write the shading parameters on every cell instance of the netlist."""
    line_re = _cell_line_regex(module)
    out, seen = [], set()
    for line in text.splitlines():
        m = line_re.match(line.strip())
        if not (m and m.group("ref") in cell_fcov):
            out.append(line)
            continue
        ref = m.group("ref")
        f = cell_fcov[ref]
        params = dict(_PARAM_TOKEN.findall(m.group("rest")))   # keep the other parameters
        for key in ("f_cov", "tau_I", "rp"):
            params.pop(key, None)

        if module.shading_params == "fcov_tau":
            new = [f"f_cov={f:.4f}", f"tau_I={tau_i}"]
        elif module.shading_params == "fcov_only":
            new = [f"f_cov={f * (1.0 - tau_i):.4f}"]
        else:
            raise ValueError(f"unknown shading_params {module.shading_params!r}")
        if module.rp_law is not None:
            # Match the coverage precision actually written into the photocurrent expression.
            new.append(f"rp={shaded_rp(module, round(f, 4), tau_i):.12g}")
        new += [f"{k}={v}" for k, v in params.items()]

        out.append(f"{ref} {m.group('n1')} {m.group('n2')} {m.group('subckt')} params: {' '.join(new)}")
        seen.add(ref)

    missing = set(cell_fcov) - seen
    if missing:
        raise RuntimeError(f"{module.name}: cells not found in {module.netlist}: {sorted(missing)[:10]}...")
    return "\n".join(out) + "\n"


def set_param(text, name, value):
    """Replace the value of `name` on the .param line(s) of the netlist."""
    token = re.compile(rf"\b({re.escape(name)})\s*=\s*([-+0-9.eE]+)", re.IGNORECASE)
    out, found = [], False
    for line in text.splitlines():
        if line.strip().lower().startswith(".param"):
            line, n = token.subn(lambda m: f"{m.group(1)}={value}", line)
            found = found or n > 0
        out.append(line)
    if not found:
        raise RuntimeError(f"no '.param {name}=...' line in the netlist")
    return "\n".join(out) + "\n"


def set_temperature(text, module, t):
    """Module temperature: `.param T=...` (cell models that use T) or the
    `.temp` directive (cell models that rely on the SPICE temperature)."""
    if module.temperature == "param":
        return set_param(text, config.TEMPERATURE_PARAM, t)
    out, found = [], False
    for line in text.splitlines():
        if line.strip().lower().startswith(".temp "):
            line, found = f".temp {t}", True
        out.append(line)
    if not found:
        raise RuntimeError("no '.temp' line in the netlist")
    return "\n".join(out) + "\n"


# --- running LTspice -------------------------------------------------------------

def run_ltspice(net_path):
    exe = config.LTSPICE_EXE
    if not exe.exists():
        raise FileNotFoundError(f"LTspice not found at {exe}")
    result = subprocess.run([str(exe), "-b", str(net_path)], capture_output=True, text=True,
                            timeout=config.LTSPICE_TIMEOUT_S)
    raw_path = net_path.with_suffix(".raw")
    if result.returncode != 0 or not raw_path.exists():
        raise RuntimeError(f"LTspice failed on {net_path.name} (code {result.returncode})\n"
                           f"{result.stdout}\n{result.stderr}")
    return raw_path


def read_iv(raw_path):
    """Module I-V curve from a .raw file, sorted by voltage, first quadrant only."""
    from PyLTSpice import RawRead          # imported here so the rest works without it
    raw = RawRead(str(raw_path))
    v = np.asarray(raw.get_trace(config.VOLTAGE_TRACE).get_wave(), dtype=float)
    i = np.asarray(raw.get_trace(config.CURRENT_TRACE).get_wave(), dtype=float) * config.CURRENT_SIGN
    order = np.argsort(v, kind="stable")
    v, i = v[order], i[order]
    keep = (v >= 0) & (i >= 0)
    return v[keep], i[keep]


def iv_parameters(v, i):
    """Isc, Voc, MPP and FF of a curve sorted by increasing voltage."""
    p = v * i
    k = int(np.argmax(p))
    isc = float(np.interp(0.0, v, i))
    voc = float(np.interp(0.0, i[::-1], v[::-1])) if i[-1] <= i[0] else float(v[-1])
    ff = 100.0 * p[k] / (isc * voc) if isc * voc else float("nan")
    return {"Isc_A": isc, "Voc_V": voc, "Imp_A": float(i[k]), "Vmp_V": float(v[k]),
            "Pmax_W": float(p[k]), "FF_percent": ff}


# --- one module, many simulations --------------------------------------------------

@dataclass
class IVResult:
    params: dict
    v: np.ndarray
    i: np.ndarray


class Simulator:
    """Runs one module. The template netlist is read once, and the clean
    reference (no band) is simulated once and cached."""

    def __init__(self, module, runs_dir=None):
        self.module = module
        self.runs_dir = runs_dir or config.RUNS_DIR
        self.runs_dir.mkdir(parents=True, exist_ok=True)
        if not module.netlist_path.exists():
            raise FileNotFoundError(f"netlist of {module.name} not found: {module.netlist_path}")
        self.template = read_netlist(module.netlist_path)
        self._done = {}            # tag -> IVResult, so a state is never simulated twice

    @classmethod
    def for_module(cls, module):
        """One simulator per module for the whole run (shared cache)."""
        key = (module.name, module.rp_law)
        if key not in _SIMULATORS:
            _SIMULATORS[key] = cls(module)
        return _SIMULATORS[key]

    def run(self, cell_fcov, tau_i, tag, temperature=None):
        """Simulate one shading state. temperature defaults to module.t_ref."""
        t = self.module.t_ref if temperature is None else temperature
        if (tag, t) in self._done:
            return self._done[(tag, t)]
        text = patch_cells(self.template, self.module, {k: round(f, 4) for k, f in cell_fcov.items()}, tau_i)
        text = set_temperature(text, self.module, t)
        net_path = self.runs_dir / f"{tag}.net"
        net_path.write_text(text, encoding="utf-8")
        v, i = read_iv(run_ltspice(net_path))
        if len(v) == 0:
            raise RuntimeError(f"empty I-V curve for {tag}")
        self._done[(tag, t)] = IVResult(iv_parameters(v, i), v, i)
        return self._done[(tag, t)]

    def clean(self):
        zero = {name: 0.0 for name in self.module.instances.values()}
        return self.run(zero, 1.0, f"clean_{self.module.name}")


_SIMULATORS = {}


def run_tag(module, shape, orientation, acov_value, tau_i, suffix=""):
    """File name of a simulation, e.g. rectangle_Acov0.06_portrait_tauI0.5_LONGI."""
    tag = f"{shape}_Acov{round(acov_value, 4)}_{orientation}_tauI{tau_i}_{module.name}"
    return f"{tag}_{suffix}" if suffix else tag
