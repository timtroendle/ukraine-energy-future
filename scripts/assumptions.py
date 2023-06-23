import pandas as pd


def assumption_table(costs: pd.DataFrame, technologies: list[str], parameters: list[str],
                     biomass_params: dict[str, str], parameter_to_format_string: dict[str, str]) -> pd.DataFrame:
    for parameter, value in biomass_params.items():
        try:
            costs.loc[("biomass", parameter), "value"] = value
        except KeyError:
            pass # This means the parameter is not in the costs df and that's fine.
    return (
        costs
        .to_xarray()
        .sel(technology=technologies, parameter=parameters)
        .value
        .to_dataframe()
        .unstack("parameter")["value"]
        .assign(**{
            col_name: format_float_col(col_name, fmt_str)
            for col_name, fmt_str in parameter_to_format_string.items()
        })
    )


def format_float_col(col: str, fmt_str: str):
    def format_col(df: pd.DataFrame):
        return df[col].map(lambda cell: f"{cell:{fmt_str}}")
    return format_col


if __name__ == "__main__":
    table = assumption_table(
        costs=pd.read_csv(snakemake.input.cost, index_col=[0, 1]),
        technologies=snakemake.params.technologies,
        parameters=snakemake.params.parameters,
        parameter_to_format_string=snakemake.params.parameter_to_format_string,
        biomass_params=snakemake.params.biomass_parameters
    )
    table.to_csv(snakemake.output[0], index=True, header=True)