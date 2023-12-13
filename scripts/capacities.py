from dataclasses import dataclass, field

import pandas as pd
import pypsa


MW_TO_GW = 1e-3

@dataclass
class Scenario:
    n: pypsa.Network
    opt_in_name: bool
    name: str = field(init=False)

    def __post_init__(self):
        self.name = self.n._meta["run"]["name"]
        if self.opt_in_name:
            self.name = self.name + "--" + self.n._meta["wildcards"]["opts"]


def all_aggregated_capacities(scenarios: list[Scenario], aggregator: callable) -> pd.DataFrame:
    return pd.DataFrame({
        s.name: aggregator(s.n)
        for s in scenarios
    })


def aggregated_power_capacities(n: pypsa.Network) -> pd.Series:
    return pd.concat([
        n.generators.groupby("carrier").p_nom_opt.sum().mul(MW_TO_GW),
        n.storage_units.groupby("carrier").p_nom_opt.sum().mul(MW_TO_GW),
        n.links.groupby("carrier").p_nom_opt.sum().mul(MW_TO_GW),
        pd.Series(index=["average-load"], data=n.loads_t.p.sum(axis=1).mul(MW_TO_GW).mean())
    ]).rename("Capacities (GW)")


def aggregated_energy_capacities(n: pypsa.Network) -> pd.Series:
    return (
        n
        .stores
        .groupby("carrier")
        .e_nom_opt
        .sum()
        .mul(MW_TO_GW)
        .rename("Capacities (GWh)")
    )


if __name__ == "__main__":
    scenarios=(
        Scenario(n=pypsa.Network(path_to_n), opt_in_name=snakemake.params.opts_out)
        for path_to_n in snakemake.input.scenarios
    )
    power = all_aggregated_capacities(scenarios, aggregated_power_capacities)
    scenarios=(
        Scenario(n=pypsa.Network(path_to_n), opt_in_name=snakemake.params.opts_out)
        for path_to_n in snakemake.input.scenarios
    )
    energy = all_aggregated_capacities(scenarios, aggregated_energy_capacities)
    power.to_csv(snakemake.output.power, index=True, header=True, float_format="%.1f")
    energy.to_csv(snakemake.output.energy, index=True, header=True, float_format="%.1f")
