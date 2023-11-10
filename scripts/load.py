import pandas as pd
import xarray as xr


def read_profile(ds: xr.Dataset, profile_name: str, scaled: bool, scale_to: str) -> pd.DataFrame:
    df = ds.to_dataframe()[snakemake.params.profile].rename("UA")
    if scaled:
        scaling_factor = calculate_scaling_factor(ds, profile_name, scale_to)
        df = df.mul(scaling_factor)
    return df


def calculate_scaling_factor(ds: xr.Dataset, profile_name: str, scale_to: str) -> float:
    # Linear scaling factor.
    total_demand = ds.sum("time")
    return (total_demand[scale_to] / total_demand[profile_name]).item()


if __name__ == "__main__":
    ds = xr.open_dataset(snakemake.input[0])
    df = read_profile(
        ds,
        profile_name=snakemake.params.profile,
        scaled=snakemake.params.scale_to is not None,
        scale_to=snakemake.params.scale_to
    )
    df.to_csv(snakemake.output[0])
