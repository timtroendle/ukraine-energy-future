import pandas as pd
import pypsa


def lcoe_all(scenarios: list[pypsa.Network], ignore_existing: bool) -> pd.Series:
    demand = demand_all(scenarios)
    total_cost = total_cost_all(scenarios, ignore_existing)
    return (total_cost / demand).rename("LCOE (€/MWh)")


def demand_all(scenarios: list[pypsa.Network]) -> pd.Series:
    return pd.Series(
        index=[s._meta["run"]["name"] for s in scenarios],
        data=[
            n.statistics.withdrawal(aggregate_time="sum").xs(("Load", "-")) * (-1)
            for n in scenarios
        ],
        name="demand"

    )


def total_cost_all(scenarios: list[pypsa.Network], ignore_existing: bool) -> pd.Series:
    return pd.Series(
        index=[s._meta["run"]["name"] for s in scenarios],
        data=[total_cost_network(s, ignore_existing) for s in scenarios],
        name="levelised cost of electricity"

    )


def total_cost_network(n: pypsa.Network, ignore_existing: bool = False) -> float:
    return sum([
        cost_generators(n, ignore_existing=ignore_existing),
        cost_storage_units(n, ignore_existing=ignore_existing),
        cost_stores(n, ignore_existing=ignore_existing),
        cost_lines(n, ignore_existing=ignore_existing),
        cost_links(n, ignore_existing=ignore_existing),
    ])


def cost_component(n: pypsa.Network, component: str, scale_factor: str, ignore_existing: bool = False) -> float:
    components = n.__getattribute__(component)
    total_cost = (components.capital_cost * components[f"{scale_factor}_opt"]).sum()
    if ignore_existing:
        existing_cost = (components.capital_cost * components[f"{scale_factor}_min"]).sum()
        return total_cost - existing_cost
    else:
        return total_cost


def cost_generators(n: pypsa.Network, ignore_existing: bool = False) -> float:
    return cost_component(n=n, component="generators", scale_factor="p_nom", ignore_existing=ignore_existing)


def cost_storage_units(n: pypsa.Network, ignore_existing: bool = False) -> float:
    return cost_component(n=n, component="storage_units", scale_factor="p_nom", ignore_existing=ignore_existing)


def cost_stores(n: pypsa.Network, ignore_existing: bool = False) -> float:
    return cost_component(n=n, component="stores", scale_factor="e_nom", ignore_existing=ignore_existing)


def cost_lines(n: pypsa.Network, ignore_existing: bool = False) -> float:
    return cost_component(n=n, component="lines", scale_factor="s_nom", ignore_existing=ignore_existing)


def cost_links(n: pypsa.Network, ignore_existing: bool = False) -> float:
    return cost_component(n=n, component="links", scale_factor="p_nom", ignore_existing=ignore_existing)


if __name__ == "__main__":
    lcoes = lcoe_all(
        scenarios=[pypsa.Network(path_to_n) for path_to_n in snakemake.input.scenarios],
        ignore_existing=snakemake.params.ignore_existing
    )
    lcoes.to_csv(snakemake.output[0], index=True, header=True, float_format="%0.1f")