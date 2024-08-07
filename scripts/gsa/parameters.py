from collections.abc import Iterable

import pandas as pd


def parameters_table(parameters: dict[str: dict[str: float]]) -> pd.DataFrame:
    return pd.DataFrame(
        index=nice_names(parameters.keys()),
        data={
            "Minimum": [f"{p['min'] * 100}%" for p in parameters.values()],
            "Maximum": [f"{p['max'] * 100}%" for p in parameters.values()]
        }
    )


def nice_names(names: Iterable[str]) -> pd.Series:
    return (
        pd
        .Series(names)
        .str
        .replace("onwind", "wind")
        .str
        .replace("biomassinflow", "Biomass")
        .str
        .replace("\+c", " capital cost", regex=True)
        .str
        .replace("\+m", " marginal cost", regex=True)
        .str
        .replace("\+l", " potential", regex=True)
        .str
        .capitalize()
    )


if __name__ == "__main__":
    table = parameters_table(
        parameters=snakemake.params.parameters
    )
    table.to_csv(snakemake.output[0], index=True, header=True)
