import pandas as pd
from SALib.sample.morris import sample as morris_sample


def input_sample(parameters: dict[str: dict[str: float]], number_trajectories: int, seed: int):
    problem = create_problem(parameters)
    param_values = morris_sample(problem, N=number_trajectories, seed=seed)
    return pd.DataFrame(param_values, columns=parameters.keys())


def create_problem(parameters: dict[str: dict[str: float]]) -> dict:
    return {
        'num_vars': len(parameters.keys()),
        'names': parameters.keys(),
        'bounds': [(param["min"], param["max"]) for param in parameters.values()]
    }


if __name__ == "__main__":
    x = input_sample(
        parameters=snakemake.params.parameters,
        number_trajectories=snakemake.params.number_trajectories,
        seed=snakemake.params.seed
    )
    x.to_csv(snakemake.output[0])
