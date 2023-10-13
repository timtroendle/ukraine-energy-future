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
SCENARIO_NAMES = [SCENARIO_NAME_1, SCENARIO_NAME_2, SCENARIO_NAME_3, SCENARIO_NAME_4]
RESULT_PATH = "pypsa-eur/results/networks/{scenario}/elec_s_{spatial_res}_ec_lvopt_{temporal_res}-BAU{opts}.nc"


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


use rule * from pypsa_scenario1 exclude build_load_data, solve_network as pypsa1_*
use rule * from pypsa_scenario2 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, solve_network as pypsa2_*
use rule * from pypsa_scenario3 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, solve_network as pypsa3_*
use rule * from pypsa_scenario4 exclude retrieve_databundle, build_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, solve_network as pypsa4_*
use rule solve_network from pypsa_scenario1 as pypsa1_solve_network with:
    resources:
        runtime = 240
use rule solve_network from pypsa_scenario2 as pypsa2_solve_network with:
    resources:
        runtime = 240
use rule solve_network from pypsa_scenario3 as pypsa3_solve_network with:
    resources:
        runtime = 240
use rule solve_network from pypsa_scenario4 as pypsa4_solve_network with:
    resources:
        runtime = 240

localrules: pypsa1_retrieve_cost_data, pypsa1_retrieve_databundle, pypsa1_retrieve_ship_raster
localrules: pypsa1_retrieve_natura_raster, pypsa1_download_copernicus_land_cover, pypsa1_build_powerplants
localrules: pypsa2_retrieve_cost_data, pypsa2_retrieve_natura_raster, pypsa2_build_powerplants
localrules: pypsa3_retrieve_cost_data, pypsa3_retrieve_natura_raster, pypsa3_build_powerplants
localrules: pypsa4_retrieve_cost_data, pypsa4_retrieve_natura_raster, pypsa4_build_powerplants


rule build_load_data_high:
    message: "Use externally defined load time series in all scenarios with high growth."
    input: "data/electricity-demand-fully-electrified-high.csv"
    output: "pypsa-eur/resources/{scenario}/load.csv"
    wildcard_constraints:
        scenario = "((nuclear-and-renewables-high)|(only-renewables-high))"
    run:
        copyfile(input[0], output[0])


rule build_load_data_low:
    message: "Use externally defined load time series in all scenarios with low growth."
    input: "data/electricity-demand-fully-electrified-low.csv"
    output: "pypsa-eur/resources/{scenario}/load.csv"
    wildcard_constraints:
        scenario = "((nuclear-and-renewables-low)|(only-renewables-low))"
    run:
        copyfile(input[0], output[0])


####################
# The following moves results from pypsa-eur to build folder of this project
####################


rule run_scenarios:
    input:
        s1 = RESULT_PATH.format(
            scenario=SCENARIO_NAME_1,
            spatial_res=config["resolution"]["space"],
            temporal_res=config["resolution"]["time"],
            opts=""),
        s2 = RESULT_PATH.format(
            scenario=SCENARIO_NAME_2,
            spatial_res=config["resolution"]["space"],
            temporal_res=config["resolution"]["time"],
            opts=""),
        s3 = RESULT_PATH.format(
            scenario=SCENARIO_NAME_3,
            spatial_res=config["resolution"]["space"],
            temporal_res=config["resolution"]["time"],
            opts=""),
        s4 = RESULT_PATH.format(
            scenario=SCENARIO_NAME_4,
            spatial_res=config["resolution"]["space"],
            temporal_res=config["resolution"]["time"],
            opts=""),
    output:
        s1 = f"build/results/scenarios/{SCENARIO_NAME_1}.nc",
        s2 = f"build/results/scenarios/{SCENARIO_NAME_2}.nc",
        s3 = f"build/results/scenarios/{SCENARIO_NAME_3}.nc",
        s4 = f"build/results/scenarios/{SCENARIO_NAME_4}.nc",
    run:
        copyfile(input.s1, output.s1)
        copyfile(input.s2, output.s2)
        copyfile(input.s3, output.s3)
        copyfile(input.s4, output.s4)


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
    x = pd.read_csv(checkpoints.gsa_input.get().output[0], index_col=0)
    pypsa_opts = {
        run_id: "-".join(f"{param_name}{x.loc[run_id, param_name]:.2f}" for param_name in x.keys())
        for run_id in x.index
    }
    runs = [
        (
            RESULT_PATH.format(
                scenario=SCENARIO_NAME_1,
                spatial_res=config["global-sensitivity-analysis"]["resolution"]["space"],
                temporal_res=config["global-sensitivity-analysis"]["resolution"]["time"],
                opts="-" + opts),
            RESULT_PATH.format(
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
        xy = "build/results/gsa/sensitivities.csv"
    conda: "../envs/gsa.yaml"
    script: "../scripts/gsa/output.py"


rule gsa_vis:
    message: "Visualise sensitivities."
    input:
        sensitivities = rules.gsa_output.output[0]
    output: "build/results/gsa/sensitivities.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/morris.py"


rule annotate_vis:
    message: "Add annotations to GSA vis."
    input:
        image = "build/results/gsa/sensitivities.png"
    output:
        image = "build/results/gsa/sensitivities-annotated.png"
    conda: "../envs/pillow.yaml"
    script: "../scripts/vis/annotate_morris.py"
