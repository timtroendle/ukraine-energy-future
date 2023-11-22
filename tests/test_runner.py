import sys
from pathlib import Path

import pytest
import pypsa


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

        @pytest.fixture(scope="session", params=snakemake.input.main_scenarios + snakemake.input.gsa_scenarios)
        def scenario(self, request):
            return pypsa.Network(request.param)

        @pytest.fixture(params=snakemake.input.main_scenario_logs + snakemake.input.gsa_scenario_logs)
        def scenario_log(self, request):
            with Path(request.param).open("r") as f_log:
                return f_log.readlines()

        @pytest.fixture(params=["offwind", "onwind", "solar"])
        def re_potential_time_series(self, scenario, request):
            all_potentials = scenario.generators_t.p_max_pu
            generators = [g for g in all_potentials.columns if request.param in g]
            return all_potentials[generators].mean(axis="columns")

        @pytest.fixture()
        def biomass_generation(self, scenario):
            periods = scenario.snapshot_weightings.generators.sum()
            return (
                scenario
                .links_t
                .p1
                .groupby(scenario.links.carrier, axis=1)
                .sum() # across locations
                .mean() # across time
                .mul(periods)
                .mul(-1)
                .loc["biomass"]
            )

        @pytest.fixture()
        def biomass_potential(self):
            return snakemake.params.biomass_potential

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(snakemake)
