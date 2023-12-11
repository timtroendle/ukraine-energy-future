import numpy as np
import pandas as pd
import xarray as xr
from SALib.analyze.morris import analyze as morris_analyze


def determine_sensitivities(parameters: dict[str: dict[str: float]], path_to_lcoe: str,
                            x: pd.DataFrame, seed: int) -> (pd.DataFrame, pd.Series):
    problem = create_problem(parameters)
    lcoes = read_lcoe(path_to_lcoe)
    Y = y(x, lcoes)
    sensitivities = (
        pd
        .DataFrame(morris_analyze(problem, x.values, Y, seed))
        .set_index("names")
        .pipe(scale_sensitivities, parameters)
    )
    return (sensitivities, pd.Series(index=x.index, data=Y, name="lcoe_diff"))


def y(x: pd.DataFrame, lcoe: xr.DataArray) -> list:
    pypsa_opts = {
        run_id: "-".join(f"{param_name}{x.loc[run_id, param_name]:.2f}" for param_name in x.keys())
        for run_id in x.index
    }
    return np.asarray([lcoe_diff(lcoe.sel(opts=opts)) for opts in pypsa_opts.values()])


def read_lcoe(path_to_file: str) -> xr.DataArray:
    return (
        pd
        .read_csv(path_to_file, index_col=0)
        .assign(
            scenario=lambda df: df.index.str.split("--").str[0],
            opts=lambda df: df.index.str.split("--").str[1].map(gsa_opts)
        )
        .drop_duplicates(["scenario", "opts"])
        .set_index(["scenario", "opts"])
        .to_xarray()
        ["LCOE (€/MWh)"]
    )


def gsa_opts(opts: str):
    opts = opts.split("-")
    assert opts[0][-1] == "H" # first option is time resolution

    first_gsa_option = 1 if opts[1] != "BAU" else 2
    return "-".join(opts[first_gsa_option:])


def lcoe_diff(lcoes: xr.DataArray) -> float:
    return (lcoes.sel(scenario="nuclear-and-renewables-high") - lcoes.sel(scenario="only-renewables-high")).item()


def create_problem(parameters: dict[str: dict[str: float]]) -> dict:
    return {
        'num_vars': len(parameters.keys()),
        'names': parameters.keys(),
        'bounds': [(param["min"], param["max"]) for param in parameters.values()]
    }


def scale_sensitivities(sensitivities: pd.DataFrame, parameters: dict[str: dict[str: float]]) -> pd.DataFrame:
    # Sensitivities are not scaled by default. We would like to report sensitivity scaled
    # to a unit change of the parameter. Thus, we need to scale be the range of the parameter.
    param_range = np.array([param["max"] - param["min"] for param in parameters.values()])
    return (sensitivities.T / param_range).T


if __name__ == "__main__":
    x = pd.read_csv(snakemake.input.x, index_col=0)
    sensitivities, lcoe_diffs = determine_sensitivities(
        parameters=snakemake.params.parameters,
        seed=snakemake.params.seed,
        x=x,
        path_to_lcoe=snakemake.input.lcoe
    )
    sensitivities.to_csv(snakemake.output.sensitivities, index=True, header=True)
    lcoe_diffs.to_csv(snakemake.output.lcoe_diffs, index=True, header=True)
    xy = pd.concat([x, lcoe_diffs], axis=1)
    xy.to_csv(snakemake.output.xy, index=True, header=True)
