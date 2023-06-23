import pandas as pd
import pypsa


MW_TO_TW = 1e-6


def all_generation_per_carrier(scenarios: list[pypsa.Network], aggregator: callable) -> pd.DataFrame:
    return pd.DataFrame({
        s._meta["run"]["name"]: aggregator(s)
        for s in scenarios
    })


def generation_per_carrier(n: pypsa.Network) -> pd.Series:
    return (
        n
        .generators_t
        .p
        .groupby(n.generators.carrier, axis=1)
        .sum() # across locations
        .mean() # across time
        .drop(index="load")
        .mul(8760)
        .mul(MW_TO_TW)
    )


if __name__ == "__main__":
    scenarios = [pypsa.Network(s) for s in snakemake.input.scenarios]
    energy = all_generation_per_carrier(scenarios, generation_per_carrier)
    energy.to_csv(snakemake.output[0], index=True, header=True, float_format="%.1f")
