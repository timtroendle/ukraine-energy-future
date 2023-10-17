import xarray as xr


if __name__ == "__main__":
    ds = xr.open_dataset(snakemake.input[0])
    df = ds.to_dataframe()[snakemake.params.profile].rename("UA")
    df.to_csv(snakemake.output[0])
