"""Discounted robot ownership costs; independent of electrical/dust models.

Currency and area conversion preserve the original investment range. Cash
outlays occur at years 0, 12, 24; maintenance is paid at each year end 1..25.
Residual treatment must be supplied explicitly by every caller.
"""
import numpy as np
import pandas as pd

EVALUATION_YEARS = 25
DISCOUNT_RATE = 0.05
ROBOT_LIFETIME_YEARS = 12
ROBOT_MAINTENANCE_FRACTION = 0.05
ROBOT_INVESTMENT_EUR_KWP = (2.4 * 4.6, 8.2 * 4.6)
RESIDUAL_MODES = ("none", "linear_remaining_life")


def annuity_factor(years=EVALUATION_YEARS, rate=DISCOUNT_RATE):
    """PV of one unit paid at the END of every year, including the final year."""
    if years <= 0 or rate <= -1:
        raise ValueError("Positive horizon and discount rate > -1 required")
    return float(np.sum((1 + rate) ** -np.arange(1, years + 1, dtype=float)))


def robot_cashflows(investment, *, residual_mode, years=EVALUATION_YEARS,
                    rate=DISCOUNT_RATE, lifetime=ROBOT_LIFETIME_YEARS,
                    maintenance_fraction=ROBOT_MAINTENANCE_FRACTION):
    """Positive net cost; negative terminal credit only in the named scenario.

    Linear residual is an explicit remaining-service-value sensitivity, NOT a
    guaranteed resale price: the final robot retains 11/12 of its original
    cost at year 25. No sale/disposal costs or escalation are assumed.
    """
    if residual_mode not in RESIDUAL_MODES:
        raise ValueError(f"residual_mode must be one of {RESIDUAL_MODES}")
    if investment < 0 or lifetime <= 0 or maintenance_fraction < 0:
        raise ValueError("Investment/maintenance must be nonnegative and life positive")
    annuity_factor(years, rate)
    t = np.arange(years + 1)
    purchases = list(range(0, years, lifetime))
    capital = np.where(np.isin(t, purchases), investment, 0.)
    maintenance = np.where(t > 0, maintenance_fraction * investment, 0.)
    credit = np.zeros_like(t, dtype=float)
    remaining_fraction = max(0., 1 - (years - purchases[-1]) / lifetime)
    if residual_mode == "linear_remaining_life":
        credit[-1] = investment * remaining_fraction
    discount = (1 + rate) ** -t.astype(float)
    cost = capital + maintenance - credit
    return pd.DataFrame(dict(year=t, investment_eur_kwp=capital,
        maintenance_eur_kwp=maintenance, residual_credit_eur_kwp=credit,
        discount_factor=discount, net_cost_eur_kwp=cost,
        discounted_net_cost_eur_kwp=cost * discount))


def robot_ownership_cost(investment, *, residual_mode):
    flows = robot_cashflows(investment, residual_mode=residual_mode)
    discount = flows.discount_factor
    pv = float(flows.discounted_net_cost_eur_kwp.sum())
    return dict(investment_eur_kwp=float(investment), residual_mode=residual_mode,
        capital_pv_eur_kwp=float((flows.investment_eur_kwp * discount).sum()),
        maintenance_pv_eur_kwp=float((flows.maintenance_eur_kwp * discount).sum()),
        residual_credit_pv_eur_kwp=float((flows.residual_credit_eur_kwp * discount).sum()),
        ownership_pv_eur_kwp=pv,
        ownership_eac_eur_kwp_year=pv / annuity_factor())


def robot_strategy_costs(rows, *, residual_mode):
    """Per-configuration totals in % of CLEAN annual revenue.

    Saved robot cost and cost_u contain energy losses only. Add ownership
    EAC exactly once in both counterfactuals. Never add per-clean manual fees.
    Aggregate the resulting ratios, not costs divided by a site's first row.
    """
    if rows.empty or not (rows["mode"] == "robot").all():
        raise ValueError("Expected nonempty robot-only rows")
    if (rows.revenue <= 0).any():
        raise ValueError("Clean annual revenue must be positive")
    # Guard the saved-input contract before adding robot ownership costs.
    np.testing.assert_allclose(rows.cost, rows.revenue * rows.loss, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(rows.cost_u, rows.revenue * rows.loss_u, rtol=1e-10, atol=1e-10)
    result = rows.copy()
    result["robot_residual_mode"] = residual_mode
    low, high = ROBOT_INVESTMENT_EUR_KWP
    for bound, initial in (("low", low), ("mid", (low + high) / 2), ("high", high)):
        eac = robot_ownership_cost(initial, residual_mode=residual_mode)["ownership_eac_eur_kwp_year"]
        result[f"ownership_eac_{bound}_eur_kwp_year"] = eac
        for case, column in (("band", "cost"), ("noband", "cost_u")):
            annual = result[column] + eac
            result[f"{case}_total_{bound}_eur_kwp_year"] = annual
            result[f"{case}_total_{bound}_pv_eur_kwp"] = annual * annuity_factor()
            result[f"{case}_{bound}_pct_clean_revenue"] = 100 * annual / result.revenue
    return result
