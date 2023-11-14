import sys
from pathlib import Path

import pytest
import xarray as xr


def run_test(snakemake):
    exit_code = pytest.main(
        [
            snakemake.input.test_dir,
            f"--html={snakemake.log[0]}",
            "--self-contained-html",
            "--verbose"
        ],
        plugins=[
            _create_config_plugin(snakemake=snakemake)
        ]
    )
    if exit_code == 0:
        Path(snakemake.output[0]).touch()
    sys.exit(exit_code)


def _create_config_plugin(snakemake):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session", params=snakemake.input.scenarios)
        def scenario(self, request):
            return xr.open_dataset(request.param)

        @pytest.fixture(params=["offwind", "onwind", "solar"])
        def re_potential_time_series(self, scenario, request):
            generators = [g.item() for g in scenario.generators_t_p_i if request.param in g.item()]
            potential = (
                scenario
                .generators_t_p_max_pu
                .sel(generators_t_p_max_pu_i=generators)
                .mean("generators_t_p_max_pu_i")
            )
            return potential

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(snakemake)
