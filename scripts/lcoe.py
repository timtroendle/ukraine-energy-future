from dataclasses import dataclass, field

import pandas as pd
import pypsa
from pypsa.descriptors import get_switchable_as_dense as as_dense


@dataclass
class Scenario:
    n: pypsa.Network
    opt_in_name: bool
    name: str = field(init=False)

    def __post_init__(self):
        self.name = self.n._meta["run"]["name"]
        if self.opt_in_name:
            self.name = self.name + "--" + self.n._meta["wildcards"]["opts"]


def lcoe_all(scenarios: list[Scenario]) -> pd.Series:
    demand = demand_all(scenarios)
    total_cost = total_cost_all(scenarios)
    return (total_cost / demand).rename("LCOE (€/MWh)")


def demand_all(scenarios: list[Scenario]) -> pd.Series:
    return pd.Series(
        index=[s.name for s in scenarios],
        data=[
            (s.n.snapshot_weightings.generators @ as_dense(s.n, "Load", "p_set")).sum()
            for s in scenarios
        ],
        name="demand"

    )


def total_cost_all(scenarios: list[Scenario]) -> pd.Series:
    return pd.Series(
        index=[s.name for s in scenarios],
        data=[total_cost_network(s.n) for s in scenarios],
        name="total system cost"

    )


def total_cost_network(n: pypsa.Network) -> float:
    return n.statistics.capex().sum() + n.statistics.opex().sum()


if __name__ == "__main__":
    lcoes = lcoe_all(
        scenarios=[Scenario(n=pypsa.Network(path_to_n), opt_in_name=snakemake.params.opts_out)
                   for path_to_n in snakemake.input.scenarios]
    )
    lcoes.to_csv(snakemake.output[0], index=True, header=True, float_format="%0.2f") # FIXME remove float format
