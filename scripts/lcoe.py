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
        .concat([lcoe(s) for s in scenarios], axis=1).T
    )


def lcoe(scenario: Scenario) -> pd.Series:
    network = scenario.n
    system = pd.Series(
        index=["System LCOE (€/MWh)"],
        data=total_cost_network(network) / demand_network(network)
    )
    components = components_lcoe(network)
    return pd.concat([system, components]).rename(scenario.name)


def components_lcoe(network: pypsa.Network) -> pd.Series():
    capex = network.statistics.capex()
    opex = network.statistics.opex().reindex(capex.index, fill_value=0)
    cost = capex + opex
    demand = demand_network(network)
    levelised_cost = cost / demand
    return levelised_cost.rename("cost").reset_index().set_index("carrier").cost


def demand_network(n: pypsa.Network) -> float:
    return (n.snapshot_weightings.generators @ as_dense(n, "Load", "p_set")).sum()


def total_cost_network(n: pypsa.Network) -> float:
    return n.statistics.capex().sum() + n.statistics.opex().sum()


def component_cost_network(n: pypsa.Network) -> pd.Series:
    return n.statistics.capex().sum() + n.statistics.opex().sum()


if __name__ == "__main__":
    lcoes = lcoe_all(
        scenarios=(Scenario(n=pypsa.Network(path_to_n), opt_in_name=snakemake.params.opts_out)
                   for path_to_n in snakemake.input.scenarios)
    )
    lcoes.to_csv(snakemake.output[0], index=True, header=True)
