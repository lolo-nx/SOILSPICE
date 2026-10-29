"""Energy yield and cost of a bottom-edge band at three sites (chapter 4).

Climate: ERA5 hourly series 2016-2025 (Open-Meteo), clean hourly DC energy from pvlib (tilt 10 deg facing the
equator, albedo 0.2, Perez transposition, physical IAM on the beam, Faiman cell temperature), with the temperature
coefficient of each module (datasheets: Yingli -0.42, SunPower -0.29, JA Solar -0.29, LONGi -0.24 %/K); the site
table and the cleaning bound use a generic module (-0.29 %/K).

Band growth (rain-event mass balance, one free parameter f):
  1. dust settles uniformly at D g m-2 per day, D = soiling rate / 4.1 % per g m-2 (Boyle et al. 2015);
  2. an event ends after three consecutive hours without rain; unless the selected criterion triggers washing,
     a total >= 1 mm washes the uniform layer down; a fraction f dries in
     the band, whose areal mass grows by f m_u L_s / h (L_s: module length along the slope);
  3. an event total >= R_H assumes complete removal at closure. The illustrative reference is R_H=15 mm/event;
     sensitivity scenarios use 5,10,20,25 mm/event. Complete washing remains uncalibrated for actual bands;
  4. tau = exp(-alpha m), alpha = 0.041 m2/g; hourly loss 1 - (1 - L_u)(1 - L_SPICE(h, tau_band)), L_SPICE from
     the equal-height sweep of chapter 2 (rectangular band), linear in tau and in h.
Cleaning strategies: none, fixed interval k days, trigger (clean when the energy lost since the last removal is
worth one cleaning), robot (daily). Prices and cleaning costs from Ilse et al. (2019).

Also (SI): washing-rain cycles with the same selected threshold, the cost of one cleaning per module
from its area, and P_an,h, the mean power while a band grows in height.
"""
from functools import lru_cache

import numpy as np
import pandas as pd

from . import config, modules, rain as rainfall
from .style import ORDER, ORIENTATIONS

# ---- sites and climate ----------------------------------------------------------------------------------------------
FILES = {"Neuchâtel": "open-meteo-46.92N7.00E519m.csv", "Almería": "open-meteo-36.87N2.49W32m.csv",
         "Doha": "open-meteo-25.34N51.43E20m.csv"}
SITES = tuple(FILES)
YEARS, TILT, ALBEDO = (2016, 2025), 10.0, 0.2
GAMMA = -0.29e-2                         # generic module (site table, cleaning bound, SI)
GAMMA_T = {"YINGLI": -0.42e-2, "SUNPOWER": -0.29e-2, "JASOLAR": -0.29e-2, "LONGI": -0.24e-2}   # datasheets, 1/K

# ---- band model -------------------------------------------------------------------------------------------------------
ALPHA = 0.041                            # m2/g, Boyle et al. 2015 (4.1 % per g m-2)
RATE = {"Neuchâtel": 0.10, "Almería": 0.28, "Doha": 0.41}   # soiling rate during dry spells, % per day
LIGHT = 1.0                             # mm/event, transfer if the event is not washing
TOTAL_THRESHOLDS = (5.0, 10.0, 15.0, 20.0, 25.0)  # mm/event, identical assumptions across scenarios
REFERENCE_TOTAL = 15.0                  # user-selected illustrative reference; not calibrated for bands
F_GRID = np.array([0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0])
F3 = [0.03, 0.1, 0.3] # f pour les calculs economiques
MARGINS = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0] 
TAU = np.linspace(0, 1, 201)             # 201 valeurs de 0 a 1 (pas de 0.005)

# ---- economics --------------------------------------------------------------------------------------------------------
MARKETS = {"ground": (0.06, 0.41), "rooftop": (0.28, 3.3)}     # EUR/kWh, EUR/kWp per cleaning (Ilse et al. 2019)
INTERVALS = [14, 30, 60, 90, 120, 180, 270, 365] # intervalles fixed de nettoyage testés
ROBOT = (2.4 * 4.6 / 10, 8.2 * 4.6 / 10)                       # min-max ; cout annuel d'un robot  EUR/kWp/yr
AF = (1 - 1.05 ** -25) / 0.05                                  # annuity factor, 25 years at 5 %; convertir un gain annuel en valeur sur 25 ans 
P_STC_W = {"YINGLI": 265.0, "SUNPOWER": 345.0, "JASOLAR": 445.0, "LONGI": 495.0}
CLEAN_EUR_M2 = {"utility": 0.09, "roof_self": 0.72}            # previous model (SI): cost of one cleaning per m2


def slope_length(m, o):
    """Module length along the slope [mm]: long side in portrait, short side in landscape."""
    return modules.MODULES[m].height if o == "portrait" else modules.MODULES[m].width


# ---- climate ------------------------------------------------------------------------------------------------------------
def era5(fname): # lit le fichier Open-Meteo 
    """Raw Open-Meteo timestamps (UTC), labelling the END of precipitation/radiation hours.
    site_hours converts to interval starts after computing the solar position at the hour midpoint.
    Source: https://open-meteo.com/en/docs/historical-weather-api#hourly-parameter-definition
    """
    path = config.CLIMATE_DIR / fname
    head = pd.read_csv(path, nrows=1)
    d = pd.read_csv(path, skiprows=3, parse_dates=["time"]).set_index("time")
    d.index = d.index.tz_localize("UTC")
    d = d.rename(columns={"precipitation (mm)": "rain", "shortwave_radiation (W/m²)": "ghi", "temperature_2m (°C)": "temp",
                          "wind_speed_10m (km/h)": "wind", "direct_normal_irradiance (W/m²)": "dni",
                          "diffuse_radiation (W/m²)": "dhi"})
    d["wind"] /= 3.6
    return d, float(head.latitude[0]), float(head.longitude[0]), float(head.elevation[0])


def pv_yield(h, lat, lon, alt, gamma=GAMMA):
    """Hourly DC energy of a clean module [kWh/kWp], temperature coefficient gamma [1/K]. Without DNI/DHI (station),
    Erbs decomposition. gamma may be a dict {name: gamma}: then one column per entry."""
    import pvlib
    sp = pvlib.solarposition.get_solarposition(h.index - pd.Timedelta("30min"), lat, lon, alt)
    sp.index = h.index
    ghi = h.ghi.clip(lower=0).fillna(0)
    if "dni" in h: # si pas DNI et DHI dans fichier .csv décomposition avec le modele
        dni, dhi = h.dni.clip(lower=0).fillna(0), h.dhi.clip(lower=0).fillna(0)
    else:
        e = pvlib.irradiance.erbs(ghi, sp.zenith, h.index)
        dni, dhi = e["dni"], e["dhi"]
    az = 180.0 if lat >= 0 else 0.0 
    poa = pvlib.irradiance.get_total_irradiance(TILT, az, sp.apparent_zenith, sp.azimuth, dni, ghi, dhi,
                                                dni_extra=pvlib.irradiance.get_extra_radiation(h.index),
                                                airmass=pvlib.atmosphere.get_relative_airmass(sp.apparent_zenith),
                                                albedo=ALBEDO, model="perez") # calcul de GPOA a 10° plein sud (Perez)
    eff = poa["poa_direct"].fillna(0) * pvlib.iam.physical(pvlib.irradiance.aoi(TILT, az, sp.apparent_zenith, sp.azimuth)) \
        + poa["poa_diffuse"].fillna(0)              # irradiance effective = direct * IAM (pertes de reflexion a angle rasant) + diffus 
    tc = pvlib.temperature.faiman(eff, h.temp.interpolate(), h.wind.interpolate().fillna(1.0)) # temperature de cellule avec le modele de Faiman (temp de l'air et vent)
    if isinstance(gamma, dict):
        return pd.DataFrame({k: (eff / 1000.0 * (1 + g * (tc - 25.0))).clip(lower=0) for k, g in gamma.items()})
    return (eff / 1000.0 * (1 + gamma * (tc - 25.0))).clip(lower=0) # return la production d'un module propre heure par heure 


@lru_cache(maxsize=None)
def site_hours(site):
    """Climate/yields on interval STARTS, retaining source years for cross-boundary rain events.
    Open-Meteo labels the preceding hour by its end. Calculate PV yield using those source labels first
    (pv_yield evaluates the sun 30 minutes earlier), then shift the joined frame once for rain/deposition.
    No boundary hours are invented when the source file does not cover them.
    """
    h, lat, lon, alt = era5(FILES[site])
    ys = pv_yield(h, lat, lon, alt, {"Y": GAMMA, **{f"Y_{m}": g for m, g in GAMMA_T.items()}})
    hours = h.join(ys)
    hours.index = hours.index - pd.Timedelta(hours=1)
    return hours


@lru_cache(maxsize=None)
def site_days(site):
    """Daily climate summaries for the site table; rain classification uses site_events instead."""
    h = site_hours(site)
    g = h.groupby(h.index.floor("D"))
    columns = [c for c in h if c == "Y" or c.startswith("Y_")]
    d = pd.DataFrame({**{c: g[c].sum() for c in columns}, "rain": g.rain.sum(),
                      "rain_max": g.rain.max()})  # largest one-hour accumulation, numerically mm/h
    return d[(d.index.year >= YEARS[0]) & (d.index.year <= YEARS[1])]


@lru_cache(maxsize=None)
def site_events(site):
    """Full event totals, including cross-midnight/cross-year rain; trailing incomplete events are retained."""
    events = rainfall.events(site_hours(site).rain)
    start = pd.Timestamp(f"{YEARS[0]}-01-01", tz="UTC")
    stop = pd.Timestamp(f"{YEARS[1] + 1}-01-01", tz="UTC")
    return events[(events.start < stop) & (events.end.isna() | (events.end > start))]


@lru_cache(maxsize=None)
def site_steps(site):
    h = site_hours(site)
    h = h[(h.index.year >= YEARS[0]) & (h.index.year <= YEARS[1])]
    return rainfall.steps(h, site_events(site))


@lru_cache(maxsize=None)
def attenuated_yields(site, D, alpha):
    """Exact sum of hourly yield * exp(-alpha D * hours/24) within each simulation interval.
    Between event closures and midnight, band mass is fixed and uniform deposition is linear in time.
    Factoring out the interval's initial uniform mass gives the same energy integral as an hourly loop.
    """
    d = site_steps(site)
    h = site_hours(site)
    h = h[(h.index.year >= YEARS[0]) & (h.index.year <= YEARS[1])]
    starts, stops = d.hour_start.to_numpy(), d.hour_stop.to_numpy()
    age_hours = np.arange(len(h)) - np.repeat(starts, stops - starts) + 1
    weights = np.exp(-alpha * D * age_hours / 24)
    columns = [c for c in h if c == "Y" or c.startswith("Y_")]
    result = np.add.reduceat(h[columns].to_numpy() * weights[:, None], starts, axis=0)
    return pd.DataFrame(result, index=d.index, columns=columns)


def washing_threshold(washing_total=REFERENCE_TOTAL):
    """Validate the assumed washing threshold in mm/event."""
    if not np.isfinite(washing_total) or washing_total <= 0:
        raise ValueError("washing threshold must be finite and positive")
    return float(washing_total)


def washing_flags(totals, *, washing_total=REFERENCE_TOTAL):
    threshold = washing_threshold(washing_total)
    values = np.asarray(totals, dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("event precipitation must be finite and nonnegative")
    return values >= threshold


def simulation_inputs(site, combos, D, alpha=ALPHA, *, washing_total=REFERENCE_TOTAL):
    """Hourly-exact inputs with a washing mask based on completed event totals."""
    washing_threshold(washing_total)
    d = site_steps(site)
    kw = dict(washing_days=washing_flags(d.rain, washing_total=washing_total),
              duration_days=d.duration_days.to_numpy(), clean_start=d.clean_start.to_numpy(),
              dust_y=yields(attenuated_yields(site, D, alpha), combos))
    return d, yields(d, combos), d.rain.to_numpy(), kw


# ---- SPICE loss of the band -------------------------------------------------------------------------------------------
@lru_cache(maxsize=None) # le calcul n'est fait que une fois 
def spice_table(shape="rectangle"): # lit equal_height_sweep.csv  il garde les 4 modules et la bande rectangulaire et prend la perte de Pmp en valeur absoule en fraction
    """Rectangular-band losses of the four modules (equal-height sweep of chapter 2): tau grid, heights and
    the loss curve against tau of every (module, orientation, height)."""
    sw = pd.read_csv(config.EQUAL_HEIGHT_CSV)
    sw = sw[sw.module.isin(ORDER) & (sw["shape"] == shape)].assign(spice=lambda x: x["Pmp_loss_%"].abs() / 100)
    ts = np.sort(sw.tau_I.unique()); hs = np.sort(sw.height_mm.unique()).astype(float) # renvoie les tau simulés  et les hauteurs simulées
    curves = {k: g.set_index("tau_I").reindex(ts).spice.to_numpy() for k, g in sw.groupby(["module", "orientation", "height_mm"])} # pour chaque (module, orientation) la perte en fonction de tau 
    return ts, hs, curves


def curve(m, o, h, shape="rectangle"): # donne la courve perte (tau) pour n'immporte quelle hauteur h de bande par interpolation lunéaire entre les deux hauteurs qui l'encadrent
    """SPICE loss against tau for a band of height h (linear in h between simulated heights, 0 at h = 0)."""
    ts, hs, curves = spice_table(shape)
    hh = np.r_[0.0, hs]; cc = np.vstack([np.zeros(len(ts)), [curves[m, o, x] for x in hs]])
    return np.array([np.interp(max(h, 0.0), hh, cc[:, i]) for i in range(len(ts))])


def lut(combos, shape="rectangle"): #construit pour chaque combinaison (module ,orientation ,h,f,g) une table de recherche de la perte sur les 201 valeurs de tau s
    """Loss look-up table in tau (TAU grid) of every combination (module, orientation, h, f, margin): with an extra
    glass margin the band keeps its height and mass, and SPICE is read at h - margin."""
    ts = spice_table(shape)[0]
    return np.array([np.interp(TAU, np.r_[ts, 1.0], np.r_[curve(m, o, h - g, shape), 0.0]) for m, o, h, f, g in combos])


def yields(d, combos):
    """Clean yield of each interval/combination, using its module's temperature coefficient."""
    return np.column_stack([d[f"Y_{c[0]}"].to_numpy() for c in combos])


def gains(combos, shape="rectangle"): # vaut f*Ls/h  (equation mb = mb + f*Ls/h)
    """Mass concentration from the glass into the band, f L_s / h."""
    k = 2.0 if shape == "triangle" else 1.0          # triangle: same trapped mass on half the area
    return np.array([k * f * slope_length(m, o) / h for m, o, h, f, g in combos])


def simulate(y, rain, D, gain, L, mode="none", k=None, price=0.06, cost=0.41, track=None,
             alpha=ALPHA, light=LIGHT, washing_days=None, *,
             washing_total=REFERENCE_TOTAL, duration_days=None, dust_y=None, clean_start=None):
    """Mass balance on intervals bounded by rain-event closures and daily cleaning decisions.
    Returns the energy lost over the period [kWh/kWp], the number of cleanings, the energy-weighted mean tau of the
    band and, with track = index, (tau, L_band, energy-weighted L_u) per interval of that combination.
    Rain contains the total of completed events at interval end; zero means no completed rainfall there.
    An explicit washing_days mask can supply the event-total decisions from simulation_inputs.
    Otherwise the mask is computed from rain and washing_total (mm/event).
    Rain acts at interval end, before any manual decision at the next interval's start (also at midnight).
    cum is a trigger bookkeeping counter, not the historical energy loss: transfer removes its uniform part,
    retaining cumb, the accumulated additional band loss. Complete washing/manual cleaning reset both.
    lost always retains all past losses. With dust_y, hourly yield and deposition are integrated exactly.
    """
    if washing_days is not None:
        washing = np.asarray(washing_days, dtype=bool)
    else:
        washing = washing_flags(rain, washing_total=washing_total)
    if washing.shape != (len(y),):
        raise ValueError("washing_days must contain one flag per simulated interval")
    dt = np.ones(len(y)) if duration_days is None else np.asarray(duration_days)
    clean_start = np.ones(len(y), bool) if clean_start is None else np.asarray(clean_start, dtype=bool)
    if dt.shape != (len(y),) or clean_start.shape != (len(y),) or not (dt > 0).all():
        raise ValueError("duration_days and clean_start must match the simulation intervals")
    if dust_y is not None and np.shape(dust_y) != np.shape(y):
        raise ValueError("dust_y must match the clean-yield array")
    n = len(gain); rr = np.arange(n)
    mu = np.zeros(n); mb = np.zeros(n); cum = np.zeros(n); cumb = np.zeros(n); age = np.zeros(n)  # mu masse uniforme et mb masse de la bande (g/m2); cum: energie perdue depuis le dernier enlevement ; cumb: partie de cum due a la bande seule
    # age : jours depuis le dernier nettoyage (pr regle interval) ; lost : energie perdue totale ; ncl: nombre de cleanings; 
    lost = np.zeros(n); ncl = np.zeros(n); tw = np.zeros(n); ser = []
    for t in range(len(y)):  # event/midnight intervals; manual decisions only at clean_start
        if not clean_start[t]:
            c = None
        elif mode == "trigger": # cleaning decisions remain once per day
            c = price * cum >= cost
        elif mode == "interval": # nettoyer si age > k jours 
            c = age >= k - 1e-9
        elif mode == "robot": # nettoyer tous les jours 
            c = np.ones(n, bool)
        else:
            c = None # ne jamais nettoyer 
        if c is not None:
            mu[c] = 0; mb[c] = 0; cum[c] = 0; cumb[c] = 0; age[c] = 0; ncl += c
        tau = np.exp(-alpha * mb)
        lb = L[rr, np.rint(tau * 200).astype(int)] # perte de la bande : lit la case de la table la plus proche de tau 
        yt = y[t]                                  # yield of the day, one value or one per combination
        yd = yt * np.exp(-alpha * D * dt[t]) if dust_y is None else dust_y[t]
        uniform_clean = np.exp(-alpha * mu) * yd
        loss = yt - (1 - lb) * uniform_clean
        mu += D * dt[t]
        lu = 1 - np.divide(uniform_clean, yt, out=np.exp(-alpha * mu).copy(), where=np.asarray(yt) > 0)
        lost += loss; cum += loss; cumb += loss - (yt - uniform_clean); tw += yt * tau; age += dt[t]
        if track is not None:
            ser.append((tau[track], lb[track], lu[track]))
        if washing[t]:  # one rain action at event closure
            mu[:] = 0; mb[:] = 0; cum[:] = 0; cumb[:] = 0; age[:] = 0
        elif rain[t] >= light:
            mb += gain * mu; mu[:] = 0; cum[:] = cumb          # the rain undid the uniform part of the cost ; la bande grandit 
    return lost, ncl, tw / y.sum(axis=0), ser


def deposition(site, factor=1.0):
    """D [g m-2 per day] from the soiling rate of the site."""
    return RATE[site] / 4.1 * factor


# ---- model runs used by the figures ------------------------------------------------------------------------------------
def model_runs(*, washing_total=REFERENCE_TOTAL):
    """Every module, orientation, h, f in F3, margin and cleaning strategy, both markets.
    Columns: site, market, mode, k, module, orientation, h, f, margin, Y, loss, loss_u, cleans, cleans_u, tau_mean,
    cost, cost_u, revenue (EUR/kWp/yr). Robot rows contain lost-revenue cost only; fig_energy adds the
    annual ROBOT investment once for its final comparison. Robot is excluded from manual-policy optima.
    """
    ts, hs, _ = spice_table()
    combos = [(m, o, h, f, g) for m in ORDER for o in ORIENTATIONS for h in hs for f in F3 for g in MARGINS]
    G, L = gains(combos), lut(combos)
    G0, L0 = np.zeros(len(combos)), np.zeros((len(combos), len(TAU)))      # uniform layer only, same module
    rows = []
    for site in SITES:
        D = deposition(site)
        d, y, rain, rain_kw = simulation_inputs(site, combos, D, washing_total=washing_total)
        nyr = d.duration_days.sum() / 365.25; Y = y.sum(axis=0) / nyr
        for mkt, (p, C) in MARKETS.items():
            runs = [("none", None), ("trigger", None)] + [("interval", k) for k in INTERVALS] + [("robot", None)]
            for mode, k in runs:
                if mkt == "rooftop" and mode not in ("none", "trigger"):
                    continue
                lb, nb, tw, _ = simulate(y, rain, D, G, L, mode, k, p, C, **rain_kw)
                l0, n0, _, _ = simulate(y, rain, D, G0, L0, mode, k, p, C, **rain_kw)
                for j, (m, o, h, f, g) in enumerate(combos):
                    rows.append(dict(site=site, market=mkt, mode=mode, k=k or 0, module=m, orientation=o, h=h, f=f,
                                     margin=g, Y=Y[j], loss=lb[j] / nyr / Y[j], loss_u=l0[j] / nyr / Y[j],
                                     cleans=nb[j] / nyr, cleans_u=n0[j] / nyr, tau_mean=tw[j]))
            print(site, mkt, "done", flush=True)
    R = pd.DataFrame(rows)
    R["cost"] = 0.0
    for mkt, (p, C) in MARKETS.items():
        k = R.market == mkt; rob = R["mode"] == "robot"
        R.loc[k, "cost"] = p * R.loss[k] * R.Y[k] + C * R.cleans[k] * (~rob[k])
        R.loc[k, "cost_u"] = p * R.loss_u[k] * R.Y[k] + C * R.cleans_u[k] * (~rob[k])
        R.loc[k, "revenue"] = p * R.Y[k]
    return R


def tau_grid(shape="rectangle", fs=F_GRID, *, washing_total=REFERENCE_TOTAL):
    """Rain only, every module/orientation/h on the full f grid: mean band tau and loss (energy maps)."""
    ts, hs, _ = spice_table()
    cf = [(m, o, h, f, 0.0) for m in ORDER for o in ORIENTATIONS for h in hs for f in fs]
    rows = []
    for site in SITES:
        d, y, rain, rain_kw = simulation_inputs(site, cf, deposition(site), washing_total=washing_total)
        nyr = d.duration_days.sum() / 365.25; Y = y.sum(axis=0) / nyr
        lb, _, tw, _ = simulate(y, rain, deposition(site), gains(cf, shape), lut(cf, shape), "none", **rain_kw)
        l0, _, _, _ = simulate(y, rain, deposition(site), np.zeros(len(cf)), np.zeros((len(cf), len(TAU))), "none", **rain_kw)
        rows += [dict(site=site, shape=shape, module=m, orientation=o, h=h, f=f, tau_mean=tw[j], loss=lb[j] / nyr / Y[j],
                      loss_u=l0[j] / nyr / Y[j]) for j, (m, o, h, f, g) in enumerate(cf)]
    return pd.DataFrame(rows)


def growth_series(ref=("JASOLAR", "landscape", 60.0, 0.1, 0.0), *, washing_total=REFERENCE_TOTAL):
    """Band state on event/midnight intervals, with energy-weighted uniform loss within each interval."""
    out = []
    for site in SITES:
        d, yy, rain, rain_kw = simulation_inputs(site, [ref], deposition(site), washing_total=washing_total)
        y = yy[:, 0]
        washing = washing_flags(d.rain, washing_total=washing_total)
        _, _, _, s = simulate(yy, rain, deposition(site), gains([ref]), lut([ref]), "none", track=0, **rain_kw)
        series = pd.DataFrame(dict(site=site, date=d.index.tz_localize(None), rain=rain, Y=y, tau=[a for a, _, _ in s],
                                   L_band=[b for _, b, _ in s], L_u=[c for _, _, c in s]))
        series["rain_max"] = d.rain_max.to_numpy()
        series["end_time"] = d.end_time.dt.tz_localize(None).to_numpy()
        series["duration_days"] = d.duration_days.to_numpy()
        series["event_id"] = d.event_id.to_numpy()
        series["washing"] = washing
        series["rain_regime"] = np.where(washing, "washing", np.where(rain >= LIGHT, "transfer",
                                          np.where(d.event_id.to_numpy() >= 0, "negligible", "no_event_end")))
        out.append(series)
    return pd.concat(out)


SENS_CASES = [("reference", {})] + [(f"{k} {s}", {k: v}) for k, lo, hi in [
    ("f", 0.03, 0.3), ("rate", 0.5, 2.0), ("light", 0.5, 2.0),
    ("alpha", 0.75, 1.25),
    ("cost", 0.5, 2.0)] for s, v in (("low", lo), ("high", hi))]


def sensitivity(h0=60.0, f0=0.1, *, washing_total=REFERENCE_TOTAL):
    """One parameter at a time around the reference (h = 60 mm, f = 0.1, ground market): energy lost to the band
    with rain only (% of yield) and extra cost of the band with the trigger rule (% of revenue), median of the 8
    configurations, and the share of modules for which landscape loses less than portrait."""
    price, cost0 = MARKETS["ground"]
    rows = []
    for case, p in SENS_CASES:
        f = f0 * (p["f"] / 0.1) if "f" in p else f0
        combos = [(m, o, h0, f, 0.0) for m in ORDER for o in ORIENTATIONS]
        G, L = gains(combos), lut(combos); G0, L0 = np.zeros(len(combos)), np.zeros_like(L)
        kw = dict(alpha=ALPHA * p.get("alpha", 1.0), light=p.get("light", LIGHT))
        cost = cost0 * p.get("cost", 1.0)
        for site in SITES:
            D = deposition(site, p.get("rate", 1.0))
            d, y, rain, rain_kw = simulation_inputs(site, combos, D, alpha=kw["alpha"], washing_total=washing_total)
            Y = y.sum(axis=0)
            lB = simulate(y, rain, D, G, L, "none", price=price, cost=cost, **rain_kw, **kw)[0]
            l0B = simulate(y, rain, D, G0, L0, "none", price=price, cost=cost, **rain_kw, **kw)[0]
            lC, nC, *_ = simulate(y, rain, D, G, L, "trigger", price=price, cost=cost, **rain_kw, **kw)
            l0C, n0C, *_ = simulate(y, rain, D, G0, L0, "trigger", price=price, cost=cost, **rain_kw, **kw)
            band = 100 * (lB - l0B) / Y
            extra = 100 * ((lC - l0C) * price + cost * (nC - n0C)) / (price * Y)
            land = np.mean([band[2 * k + 1] < band[2 * k] for k in range(len(ORDER))])
            rows.append(dict(site=site, case=case, band_loss=np.median(band), band_cost=np.median(extra),
                             landscape_better=land))
        print(case, "done", flush=True)
    return pd.DataFrame(rows)


def cleaning_bound(*, washing_total=REFERENCE_TOTAL):
    """Smallest band loss for which one cleaning can pay, L_min = (C/p) / E_next, E_next the mean energy produced
    until the next event with total precipitation >= washing_total."""
    rows = []
    for site in SITES:
        d = site_steps(site); y = d.Y.to_numpy()
        th = washing_threshold(washing_total)
        for th in (th,):
            washing = washing_flags(d.rain, washing_total=washing_total)
            En = np.full(len(y), np.nan); acc = 0.0; seen = False
            for t in range(len(y) - 1, -1, -1):
                if washing[t]:
                    acc = 0.0; seen = True
                acc += y[t]
                if seen:
                    En[t] = acc
            # Manual cleaning is considered at midnight, not at every weather-grid split.
            candidates = En[d.clean_start.to_numpy()]
            E = np.nanmean(candidates) if np.isfinite(candidates).any() else np.nan
            for mkt, (p, C) in MARKETS.items():
                rows.append(dict(site=site, washing_event_total_mm=th,
                                 washing_per_yr=washing.sum() / d.index.year.nunique(), market=mkt, E_next=E,
                                 L_min=100 * C / p / E if np.isfinite(E) else np.nan))
    return pd.DataFrame(rows)


def site_table(*, washing_total=REFERENCE_TOTAL):
    """Table sites: ERA5 grid point, yield (generic module, -0.29 %/K), rain, rainy days, soiling rate and deposition."""
    rows = []
    for site in SITES:
        d = site_days(site); ny = d.index.year.nunique()
        steps = site_steps(site)
        washing = washing_flags(steps.rain, washing_total=washing_total)
        head = pd.read_csv(config.CLIMATE_DIR / FILES[site], nrows=1)          # ERA5 grid point
        rows.append({"site": site, "Latitude": float(head.latitude[0]), "Longitude": float(head.longitude[0]),
                     "Yield [kWh/kWp/yr]": d.Y.sum() / ny, "Precipitation [mm/yr]": d.rain.sum() / ny,
                     "Days >= 1 mm [/yr]": (d.rain >= LIGHT).sum() / ny,
                     "Washing events [/yr]": washing.sum() / ny,
                     "Transfer events [/yr]": ((steps.rain >= LIGHT) & ~washing).sum() / ny,
                     "Soiling rate [%/day]": RATE[site], "D [g/m2/day]": deposition(site)})
    return pd.DataFrame(rows)


# ---- cycle-based inverse diagnosis and SI --------------------------------------------------------------------------------
def clean_cost_kwp(module_name, eur_m2):
    """Cost of one cleaning [EUR/kWp] from its cost per m2 of module."""
    m = modules.get(module_name)
    return eur_m2 * (m.width * m.height / 1e6) / (P_STC_W[module_name] / 1000.0)


def cleaning_kwh(market="utility", price=0.06):
    """Cost of one cleaning in kWh/kWp (C/p), median of the four modules."""
    return float(np.median([clean_cost_kwp(m, CLEAN_EUR_M2[market]) / price for m in ORDER]))


def washing_cycles(*, washing_total=REFERENCE_TOTAL):
    """Mean clean energy per assumed washing event (total >= R_H) [kWh/kWp]."""
    out = {}
    for site in SITES:
        d = site_steps(site)
        n_washing = washing_flags(d.rain, washing_total=washing_total).sum()
        out[site] = d.Y.sum() / n_washing if n_washing else np.inf
    return out


def pan_h(taus=np.round(np.arange(0, 0.95, 0.1), 1), hmax=np.arange(10, 101, 1.0)):
    """Mean normalised power while a rectangular band grows linearly in height from 0 to h_max."""
    sweep = pd.read_csv(config.EQUAL_HEIGHT_CSV)
    sweep = sweep[sweep["shape"] == "rectangle"]
    rows = []
    for m in ORDER:
        for o in ORIENTATIONS:
            for tau in taus:
                s = sweep[(sweep.module == m) & (sweep.orientation == o) & np.isclose(sweep.tau_I, tau)].sort_values("height_mm")
                h, L = np.r_[0, s.height_mm], np.r_[0, s["Pmp_loss_%"].abs() / 100]
                for hm in hmax:
                    pan = 1 - np.interp(np.linspace(0, hm, 401), h, L).mean()
                    rows.append(dict(module=m, orientation=o, tau_I=tau, h_max_mm=hm, P_an_h=pan,
                                     cycle_loss_pct=100 * (1 - pan), final_loss_pct=100 * np.interp(hm, h, L)))
    return pd.DataFrame(rows)


def cycle_loss(tau, loss_row, tau_f, clean="none"):
    """Loss over a cycle in which the band darkens as tau_f ** (t/T) (L interpolated in tau, L(1) = 0).
    clean: none -> cycle loss; opt -> loss saved by one cleaning at the best time (searched every T/40)."""
    s = np.linspace(0, 1, 2001); tcs = np.linspace(0.025, 0.975, 39)

    def L(x):
        return np.interp(max(tau_f, 1e-3) ** x, np.r_[tau, 1.0], np.r_[loss_row, 0.0])
    full = L(s).mean()
    if clean != "opt" or full == 0:
        return full
    ls = L(s)
    return full - min(np.where(s < tc, ls, L(np.clip(s - tc, 0, 1))).mean() for tc in tcs)


# ---- design lever 1: band height at 1 % loss ---------------------------------------------------------------------------
def frame_to_cell(m, o):
    """Distance from the frame edge to the first cell on the bottom edge [mm]."""
    M = modules.MODULES[m]
    return M.margins[2] + M.frame[2] if o == "portrait" else M.margins[0] + M.frame[0]


def design_h1(level=1.0, shapes=("rectangle", "triangle")):
    """Band height (from the frame edge) at which the SPICE loss first reaches `level` %, for every module,
    orientation and tau_I (rectangular band, linear between simulated heights; nan if never)."""
    sw = pd.read_csv(config.EQUAL_HEIGHT_CSV)
    sw = sw[sw["shape"].isin(shapes) & sw.module.isin(ORDER)].assign(L=lambda x: x["Pmp_loss_%"].abs())

    def h_at(q):
        q = q.sort_values("height_mm"); h = q.height_mm.to_numpy(); L = q.L.to_numpy()
        if not (L >= level).any():
            return np.nan
        i = int(np.argmax(L >= level))
        return np.interp(level, L[:i + 1], h[:i + 1]) if i else h[0]
    return pd.DataFrame([dict(shape=sh, module=m, orientation=o, tau_I=t, d=frame_to_cell(m, o), h1=h_at(q))
                         for (sh, m, o, t), q in sw.groupby(["shape", "module", "orientation", "tau_I"])])
