import copy
from itertools import chain
from shutil import copyfile
from pathlib import Path

from snakemake.io import load_configfile
from snakemake.utils import update_config
import pandas as pd

base_pypsa_config = load_configfile("pypsa-eur/config.default.yaml")

SCENARIO_NAME_1 = "nuclear-and-renewables-high"
SCENARIO_NAME_2 = "only-renewables-high"
SCENARIO_NAME_3 = "nuclear-and-renewables-low"
SCENARIO_NAME_4 = "only-renewables-low"
SCENARIO_NAME_5 = "only-renewables-scaled"
SCENARIO_NAME_6 = "nuclear-and-renewables-today"
SCENARIO_NAMES = [SCENARIO_NAME_1, SCENARIO_NAME_2, SCENARIO_NAME_3, SCENARIO_NAME_4, SCENARIO_NAME_5, SCENARIO_NAME_6]
RESULT_PATH = "pypsa-eur/results/networks/{scenario}/elec_s_{spatial_res}_ec_lvopt_{temporal_res}-BAU{opts}.nc"
LOG_PATH = "pypsa-eur/logs/{scenario}/solve_network/elec_s_{spatial_res}_ec_lvopt_{temporal_res}-BAU{opts}_python.log"

RUNTIME_RUNS = 600 # mins
MEMORY_RUNS = 16000 # MB

def build_scenario_config(base_pypsa, scenario_name):
    scenario_config = copy.deepcopy(base_pypsa)
    update_config(scenario_config, config["pypsa-eur"])
    update_config(scenario_config, config["scenarios"][scenario_name])
    return scenario_config


scenario_configs = [
    build_scenario_config(base_pypsa_config, scenario_name) for scenario_name in SCENARIO_NAMES
]


module pypsa_scenario1:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario_configs[0]


module pypsa_scenario2:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario_configs[1]


module pypsa_scenario3:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario_configs[2]


module pypsa_scenario4:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario_configs[3]


module pypsa_scenario5:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario_configs[4]


module pypsa_scenario6:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario_configs[5]


use rule * from pypsa_scenario1 exclude build_load_data, build_cutout, build_renewable_profiles, solve_network as pypsa1_*
use rule * from pypsa_scenario2 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, build_ship_raster, build_renewable_profiles, solve_network as pypsa2_*
use rule * from pypsa_scenario3 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, build_ship_raster, build_renewable_profiles, solve_network as pypsa3_*
use rule * from pypsa_scenario4 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, build_ship_raster, build_renewable_profiles, solve_network as pypsa4_*
use rule * from pypsa_scenario5 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, build_ship_raster, build_renewable_profiles, solve_network as pypsa5_*
use rule * from pypsa_scenario6 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, build_ship_raster, build_renewable_profiles, solve_network as pypsa6_*
use rule build_cutout from pypsa_scenario1 as pypsa1_build_cutout with: # TODO requires internet access
                                                                        # on compute notes which is not
                                                                        # given by default on Euler.
                                                                        # Load env module manually.
    resources:
        runtime = 2880
use rule solve_network from pypsa_scenario1 as pypsa1_solve_network with:
    resources:
        runtime = RUNTIME_RUNS,
        mem_mb = MEMORY_RUNS
use rule solve_network from pypsa_scenario2 as pypsa2_solve_network with:
    resources:
        runtime = RUNTIME_RUNS,
        mem_mb = MEMORY_RUNS
use rule solve_network from pypsa_scenario3 as pypsa3_solve_network with:
    resources:
        runtime = RUNTIME_RUNS,
        mem_mb = MEMORY_RUNS
use rule solve_network from pypsa_scenario4 as pypsa4_solve_network with:
    resources:
        runtime = RUNTIME_RUNS,
        mem_mb = MEMORY_RUNS
use rule solve_network from pypsa_scenario5 as pypsa5_solve_network with:
    resources:
        runtime = RUNTIME_RUNS,
        mem_mb = MEMORY_RUNS
use rule solve_network from pypsa_scenario6 as pypsa6_solve_network with:
    resources:
        runtime = RUNTIME_RUNS,
        mem_mb = MEMORY_RUNS
use rule build_renewable_profiles from pypsa_scenario1 as pypsa1_build_renewable_profiles with:
    resources:
        runtime = 60
use rule build_renewable_profiles from pypsa_scenario2 as pypsa2_build_renewable_profiles with:
    resources:
        runtime = 60
use rule build_renewable_profiles from pypsa_scenario3 as pypsa3_build_renewable_profiles with:
    resources:
        runtime = 60
use rule build_renewable_profiles from pypsa_scenario4 as pypsa4_build_renewable_profiles with:
    resources:
        runtime = 60
use rule build_renewable_profiles from pypsa_scenario5 as pypsa5_build_renewable_profiles with:
    resources:
        runtime = 60
use rule build_renewable_profiles from pypsa_scenario6 as pypsa6_build_renewable_profiles with:
    resources:
        runtime = 60


localrules: pypsa1_retrieve_cost_data, pypsa1_retrieve_databundle, pypsa1_retrieve_ship_raster
localrules: pypsa1_retrieve_natura_raster, pypsa1_download_copernicus_land_cover, pypsa1_build_powerplants
localrules: pypsa2_retrieve_cost_data, pypsa2_retrieve_natura_raster, pypsa2_build_powerplants
localrules: pypsa3_retrieve_cost_data, pypsa3_retrieve_natura_raster, pypsa3_build_powerplants
localrules: pypsa4_retrieve_cost_data, pypsa4_retrieve_natura_raster, pypsa4_build_powerplants
localrules: pypsa5_retrieve_cost_data, pypsa5_retrieve_natura_raster, pypsa5_build_powerplants
localrules: pypsa6_retrieve_cost_data, pypsa6_retrieve_natura_raster, pypsa6_build_powerplants


rule build_load_data_high:
    message: "Use externally defined load time series in all scenarios with high growth."
    input: "data/electricity-demand-fully-electrified.nc"
    params:
        profile = "2060 High Growth",
        scale_to = None
    output: "pypsa-eur/resources/{scenario}/load.csv"
    wildcard_constraints:
        scenario = "((nuclear-and-renewables-high)|(only-renewables-high)|(nuclear-and-renewables-today))"
    conda: "../envs/default.yaml"
    script: "../scripts/load.py"



rule build_load_data_low:
    message: "Use externally defined load time series in all scenarios with low growth."
    input: "data/electricity-demand-fully-electrified.nc"
    params:
        profile = "2060 Low Growth",
        scale_to = None
    output: "pypsa-eur/resources/{scenario}/load.csv"
    wildcard_constraints:
        scenario = "((nuclear-and-renewables-low)|(only-renewables-low))"
    conda: "../envs/default.yaml"
    script: "../scripts/load.py"


rule build_load_data_scaled:
    message: "Use externally defined load time series in all scenarios with scaled growth."
    input: "data/electricity-demand-fully-electrified.nc"
    params:
        profile = "2060 Low Growth",
        scale_to = "2060 High Growth"
    output: "pypsa-eur/resources/{scenario}/load.csv"
    wildcard_constraints:
        scenario = "((only-renewables-scaled))"
    conda: "../envs/default.yaml"
    script: "../scripts/load.py"


rule build_ship_raster: # this is necessary because of race conditions across scenarios
    message: "Copy existing ship raster."
    input: f"pypsa-eur/resources/{SCENARIO_NAME_1}/shipdensity_raster.nc"
    output: "pypsa-eur/resources/{scenario}/shipdensity_raster.nc"
    wildcard_constraints:
        scenario = f"(({SCENARIO_NAME_2})|({SCENARIO_NAME_3})|({SCENARIO_NAME_4})|({SCENARIO_NAME_5})|({SCENARIO_NAME_6}))"
    run:
        copyfile(input[0], output[0])


####################
# The following moves results from pypsa-eur to build folder of this project
####################


rule run_scenarios:
    input:
        scenarios = [
            RESULT_PATH.format(
                scenario=scenario,
                spatial_res=config["resolution"]["space"],
                temporal_res=config["resolution"]["time"],
                opts=""
            )
            for scenario in SCENARIO_NAMES
        ],
        logs = [
            LOG_PATH.format(
                scenario=scenario,
                spatial_res=config["resolution"]["space"],
                temporal_res=config["resolution"]["time"],
                opts=""
            )
            for scenario in SCENARIO_NAMES
        ]
    output:
        scenarios = [
            f"build/results/scenarios/{scenario}.nc"
            for scenario in SCENARIO_NAMES
        ],
        logs = [
            f"build/results/scenarios/{scenario}.log"
            for scenario in SCENARIO_NAMES
        ]
    run:
        for s_origin, s_destination in zip(input.scenarios, output.scenarios):
            copyfile(s_origin, s_destination)
        for l_origin, l_destination in zip(input.logs, output.logs):
            copyfile(l_origin, l_destination)


checkpoint gsa_input:
    message: "Create input samples for GSA."
    params:
        parameters = config["global-sensitivity-analysis"]["parameters"],
        number_trajectories = config["global-sensitivity-analysis"]["number-trajectories"],
        seed = config["global-sensitivity-analysis"]["seed"]
    output:
        x = "build/results/gsa/input.csv"
    conda: "../envs/gsa.yaml"
    script: "../scripts/gsa/input.py"


def gsa_runs(wildcards) -> list[str]:
    """Returns list with pathnames of all necessary runs for GSA."""
    return gsa_files(wildcards, RESULT_PATH)


def gsa_run_logs(wildcards) -> list[str]:
    """Returns list with pathnames of all solver logs of GSA."""
    return gsa_files(wildcards, LOG_PATH)


def gsa_files(wildcards, base_path):
    x = pd.read_csv(checkpoints.gsa_input.get().output[0], index_col=0)
    pypsa_opts = {
        run_id: "-".join(f"{param_name}{x.loc[run_id, param_name]:.2f}" for param_name in x.keys())
        for run_id in x.index
    }
    runs = [
        (
            base_path.format(
                scenario=SCENARIO_NAME_1,
                spatial_res=config["global-sensitivity-analysis"]["resolution"]["space"],
                temporal_res=config["global-sensitivity-analysis"]["resolution"]["time"],
                opts="-" + opts),
            base_path.format(
                scenario=SCENARIO_NAME_2,
                spatial_res=config["global-sensitivity-analysis"]["resolution"]["space"],
                temporal_res=config["global-sensitivity-analysis"]["resolution"]["time"],
                opts="-" + opts),
        )
        for run_id, opts in pypsa_opts.items()
    ]
    return list(chain.from_iterable(runs))


rule gsa_lcoe:
    message: "Calculate LCOE of all GSA scenarios."
    input:
        scenarios = gsa_runs
    output: "build/results/gsa/lcoe.csv"
    params:
        ignore_existing = False,
        opts_out = True
    resources:
        mem_mb = 16000
    conda: "../envs/default.yaml"
    script: "../scripts/lcoe.py"


rule gsa_output:
    message: "Analyse sensitivities of GSA."
    input:
        x = "build/results/gsa/input.csv",
        lcoe = rules.gsa_lcoe.output[0]
    params:
        parameters = config["global-sensitivity-analysis"]["parameters"],
        seed = config["global-sensitivity-analysis"]["seed"]
    output:
        sensitivities = "build/results/gsa/sensitivities.csv",
        lcoe_diffs = "build/results/gsa/output.csv",
        xy = "build/results/gsa/xy.csv"
    conda: "../envs/gsa.yaml"
    script: "../scripts/gsa/output.py"


rule gsa_vis:
    message: "Visualise sensitivities."
    input:
        sensitivities = rules.gsa_output.output.sensitivities
    output: "build/results/gsa/sensitivities.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/morris.py"


rule annotate_vis:
    message: "Add annotations to GSA vis."
    input:
        image = "build/results/gsa/sensitivities.png"
    params:
        zero_line = config["report"]["global-sensitivity-analysis"]["zero-line"]
    output:
        image = "build/results/gsa/sensitivities-annotated.png"
    conda: "../envs/pillow.yaml"
    script: "../scripts/vis/annotate_morris.py"
