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
    return (
        pd
        .DataFrame
        .from_dict(
            dict([lcoe(s) for s in scenarios]),
            orient="index",
            columns=["LCOE (€/MWh)"]
        )
    )


def lcoe(scenario: Scenario) -> (str, float):
    network = scenario.n
    return scenario.name, total_cost_network(network) / demand_network(network)


def demand_network(n: pypsa.Network) -> float:
    return (n.snapshot_weightings.generators @ as_dense(n, "Load", "p_set")).sum()


def total_cost_network(n: pypsa.Network) -> float:
    return n.statistics.capex().sum() + n.statistics.opex().sum()


if __name__ == "__main__":
    lcoes = lcoe_all(
        scenarios=(Scenario(n=pypsa.Network(path_to_n), opt_in_name=snakemake.params.opts_out)
                   for path_to_n in snakemake.input.scenarios)
    )
    lcoes.to_csv(snakemake.output[0], index=True, header=True, float_format="%0.2f") # FIXME remove float format
