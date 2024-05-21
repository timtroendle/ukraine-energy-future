import pandas as pd

VARIANT_TYPE = pd.CategoricalDtype(
    categories=[
        "Historical",
        "Low",
        "Medium",
        "High",
    ],
    ordered=True,
)

RENAMED_CATEGORIES = [
    "Historical",
    "Low fertility scenario",
    "Medium fertility scenario",
    "High fertility scenario",
]


def population(
    path_to_medium_population: str, path_to_high_low_population: str
) -> pd.DataFrame:
    medium = pd.read_csv(
        path_to_medium_population,
        low_memory=False,
        usecols=["Location", "Variant", "Time", "TPopulation1Jan"],
    )
    high_low = pd.read_csv(
        path_to_high_low_population,
        low_memory=False,
        usecols=["Location", "Variant", "Time", "TPopulation1Jan"],
    )
    all = pd.concat([
        medium,
        high_low[(high_low.Variant == "High") | (high_low.Variant == "Low")],
    ])
    ukraine = all[(all.Location == "Ukraine") & (all.Time <= 2060)].assign(
        Time=lambda df: pd.to_datetime(df.Time, format="%Y"),
        population=lambda df: df.TPopulation1Jan / 1000,
        Variant=lambda df: df.Variant.astype(VARIANT_TYPE).cat.rename_categories(
            RENAMED_CATEGORIES
        ),
    )
    historic = ukraine[
        (ukraine.Variant == "Medium fertility scenario") & (ukraine.Time <= "2022")
    ].assign(Variant="Historical")
    return pd.concat([ukraine, historic]).reset_index(drop=True)


if __name__ == "__main__":
    df_pop = population(
        path_to_medium_population=snakemake.input.medium,
        path_to_high_low_population=snakemake.input.high_low,
    )
    df_pop.to_feather(snakemake.output[0])
