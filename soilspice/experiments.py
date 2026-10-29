
"""Experimental results.

Band with height h = 80 mm applied on each module in portrait and landscape.
SUSPI transmittance values:
    1 layer  : tau_I = 0.706
    2 layers : tau_I = 0.536
    3 layers : tau_I = 0.442
    Opaque   : tau_I = 0.0
"""

from . import config

BAND_HEIGHT_MM = 80.0
TRIANGLE_APEX = "right"

_L = config.DATA_DIR / "Longi"
_S = config.DATA_DIR / "SUNPOWER"
_Y = config.DATA_DIR / "Yingli"
_J = config.DATA_DIR / "JASolar"
_S2 = config.DATA_DIR / "SPR-2"  # Second SunPower module, opaque case only


def _flashes(folder, *pairs):
    return [
        {"path": folder / name, "label": f"tau_I={tau:g}", "tau_I": tau}
        for name, tau in pairs
    ]


CAMPAIGN = {
    ("LONGI", "rectangle", "portrait"): _flashes(
        _L, ("Longi-clean-flash1-25.asc", 1.0),
        ("Longi-1L-rectangle-portrait-25.7-flash1.asc", 0.706),
        ("Longi-2L-rectangle--portrait-25.7.asc", 0.536),
        ("Longi-3L-rectangle--portrait-25.9.asc", 0.442),
        ("Longi-carton-rectangle-portrait-25.7.asc", 0.0)),

    ("LONGI", "rectangle", "landscape"): _flashes(
        _L, ("Longi-clean-flash1-25.asc", 1.0),
        ("Longi-1L-rectangle-landscape-24.7.asc", 0.706),
        ("Longi-2L-rectangle-landscape-25.4.asc", 0.536),
        ("Longi-3L-rectangle-landscape-25.9.asc", 0.442),
        ("Longi-carton-rectangle-landscape-25.5.asc", 0.0)),

    ("LONGI", "triangle", "portrait"): _flashes(
        _L, ("Longi-clean-flash1-25.asc", 1.0),
        ("Longi-1L-triangle-portrait-26.2-flash1.asc", 0.706),
        ("Longi-2L-triangle--portrait-26.3.asc", 0.536),
        ("Longi-3L-triangle--portrait-26.1.asc", 0.442),
        ("Longi-carton-triangle-portrait-25.8.asc", 0.0)),

    ("LONGI", "triangle", "landscape"): _flashes(
        _L, ("Longi-clean-flash1-25.asc", 1.0),
        ("Longi-1L-triangle-landscape-25.asc", 0.706),
        ("Longi-2L-triangle-landscape-25.4.asc", 0.536),
        ("Longi-3L-triangle-landscape-26.5.asc", 0.442),
        ("Longi-carton-triangle-landscape-25.7.asc", 0.0)),

    ("SUNPOWER", "rectangle", "portrait"): _flashes(
        _S, ("SUNPOWER-flash1-clean-25.4.asc", 1.0),
        ("SUNPOWER-flash1-26-rectangle-portrait-1layer-80m.asc", 0.706),
        ("SUNPOWER-flash1-26.2-rectangle-portrait-2layer-80m.asc", 0.536),
        ("SUNPOWER-flash1-26-rectangle-portrait-3layer-80m.asc", 0.442)) + _flashes(
        _S2, ("SPR-2-carton-portrait-rectangle-24.9.asc", 0.0)),

    ("SUNPOWER", "rectangle", "landscape"): _flashes(
        _S, ("SUNPOWER-flash1-clean-25.4.asc", 1.0),
        ("SUNPOWER-flash2-25.5-rectangle-landscape-1layer-80m.asc", 0.706),
        ("SUNPOWER-flash1-25.8-rectangle-landscape-2layer-80m.asc", 0.536)) + _flashes(
        _S2, ("SPR-2-carton-landscape-rectangle-25.3.asc", 0.0)),

    ("SUNPOWER", "triangle", "landscape"): _flashes(
        _S, ("SUNPOWER-flash1-clean-25.4.asc", 1.0),
        ("SUNPOWER-flash2-25.9-triangle-landscape-1layer-80m.asc", 0.706),
        ("SUNPOWER-flash1-25.7-triangle-landscape-2layer-80m.asc", 0.536)) + _flashes(
        _S2, ("SPR-2-carton-landscape-triangle-25.9.asc", 0.0)),

    ("SUNPOWER", "triangle", "portrait"): _flashes(
        _S, ("SUNPOWER-flash1-clean-25.4.asc", 1.0),
        ("SUNPOWER-flash1-26.2-triangle-portrait-1layer-80m.asc", 0.706),
        ("SUNPOWER-flash1-25.9-triangle-portrait-2layer-80m.asc", 0.536),
        ("SUNPOWER-flash1-26.3-triangle-portrait-3layer-80m.asc", 0.442)) + _flashes(
        _S2, ("SPR-2-carton-portrait-triangle-25.2.asc", 0.0)),

    ("YINGLI", "rectangle", "landscape"): _flashes(
        _Y, ("flash2-yingli-24.8-clean.asc", 1.0),
        ("flash1-yingli-25.9-landscape-rectangle-80m-1layer.asc", 0.706),
        ("flash1-yingli-25.4-landscape-rectangle-80m-2layer.asc", 0.536),
        ("08-09/YS-0809-carton-rectangle-landscape-24.4.asc", 0.0)),

    ("YINGLI", "rectangle", "portrait"): _flashes(
        _Y, ("flash2-yingli-24.8-clean.asc", 1.0),
        ("flash4-yingli-25.8-portrait-rectangle-80m.asc", 0.706),
        ("flash1-yingli-26.2-portrait-rectangle-80m-2layer.asc", 0.536),
        ("flash1-yingli-26.5-portrait-rectangle-80m-3layer.asc", 0.442),
        ("08-09/YS-0809-carton-rectangle-portrait-25.5.asc", 0.0)),

    ("YINGLI", "triangle", "portrait"): _flashes(
        _Y, ("flash2-yingli-24.8-clean.asc", 1.0),
        ("flahs2-yingli-26.2-triangle-portrait-80m-1layer.asc", 0.706),
        ("flahs2-yingli-26-triangle-portrait-80m-2layer.asc", 0.536),
        ("flahs1-yingli-25.8-triangle-portrait-80m-3layer.asc", 0.442),
        ("08-09/YS-0809-carton-triamgle-portrait-25.9.asc", 0.0)),

    ("YINGLI", "triangle", "landscape"): _flashes(
        _Y, ("flash2-yingli-24.8-clean.asc", 1.0),
        ("flash2-yingli-25.6-landscape-triangle-80m-1layer.asc", 0.706),
        ("flash1-yingli-25.6-landscape-triangle-80m-2layer.asc", 0.536),
        ("08-09/YS-0809-carton-triangle-landscape-25.3.asc", 0.0)),

    ("JASOLAR", "rectangle", "portrait"): _flashes(
        _J, ("JASOLAR-clean-flash1-23.7.asc", 1.0),
        ("JASOLAR-1L-rectangle-portrait-23.8.asc", 0.706),
        ("JASOLAR-2L-rectangle-portrait-24.3.asc", 0.536),
        ("JASOLAR-3L-rectangle-portrait-25.1.asc", 0.442),
        ("JASOLAR-carton-recatngle-portrait-25.1.asc", 0.0)),

    ("JASOLAR", "rectangle", "landscape"): _flashes(
        _J, ("JASOLAR-clean-flash1-23.7.asc", 1.0),
        ("JASOLAR-1L-Rectangle-Landscape-23.7.asc", 0.706),
        ("JASOLAR-2L-rectangle-Landscape-24.asc", 0.536),
        ("JASOLAR-3L-rectangle-Landscape-24.5.asc", 0.442),
        ("JASOLAR-carton-recatngle-landscape-25.asc", 0.0)),

    ("JASOLAR", "triangle", "portrait"): _flashes(
        _J, ("JASOLAR-clean-flash1-23.7.asc", 1.0),
        ("JASOLAR-1L-triangle-portrait-24.1.asc", 0.706),
        ("JASOLAR-2L-triangle-portrait-25.asc", 0.536),
        ("JASOLAR-3L-triangle-portrait-25.8.asc", 0.442),
        ("JASOLAR-carton-triangle-portrait-25.1.asc", 0.0)),

    ("JASOLAR", "triangle", "landscape"): _flashes(
        _J, ("JASOLAR-clean-flash1-23.7.asc", 1.0),
        ("JASOLAR-1L-triangle-Landscape-23.7.asc", 0.706),
        ("JASOLAR-2L-triangle-Landscape-24.3.asc", 0.536),
        ("JASOLAR-3L-triangle-Landscape-24.5.asc", 0.442),
        ("JASOLAR-carton-triangle-landscape-25.2.asc", 0.0)),
}


def flashes(module_name, shape, orientation):
    entries = CAMPAIGN.get((module_name, shape, orientation), [])
    return sorted(entries, key=lambda e: e["tau_I"], reverse=True)


def measured_modules(order):
    """Modules of `order` that have flash tests in the campaign."""
    names = {key[0] for key in CAMPAIGN}
    return [name for name in order if name in names]
