import copy
from itertools import chain
from shutil import copyfile
from pathlib import Path

from snakemake.io import load_configfile
from snakemake.utils import update_config
import pandas as pd

base_pypsa_config = load_configfile("pypsa-eur/config.default.yaml")

SCENARIO_NAME_1 = "nuclear-and-renewables"
SCENARIO_NAME_2 = "only-renewables"
RESULT_PATH = "pypsa-eur/results/networks/{scenario}/elec_s_{spatial_res}_ec_lvopt_{temporal_res}-BAU{opts}.nc"

scenario1_config = copy.deepcopy(base_pypsa_config)
scenario2_config = copy.deepcopy(base_pypsa_config)
update_config(scenario1_config, config["pypsa-eur"])
update_config(scenario1_config, config["scenarios"][SCENARIO_NAME_1])
update_config(scenario2_config, config["pypsa-eur"])
update_config(scenario2_config, config["scenarios"][SCENARIO_NAME_2])


module pypsa_scenario1:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario1_config


module pypsa_scenario2:
    snakefile:
        "../pypsa-eur/Snakefile"
    prefix:
        "pypsa-eur"
    config:
        scenario2_config


use rule * from pypsa_scenario1 exclude build_load_data as pypsa1_*
use rule * from pypsa_scenario2 exclude retrieve_databundle, retrieve_cutout, download_copernicus_land_cover, retrieve_load_data, build_load_data, retrieve_ship_raster, build_ship_raster as pypsa2_*


localrules: pypsa1_retrieve_cost_data, pypsa1_retrieve_cutout, pypsa1_retrieve_databundle, pypsa1_retrieve_ship_raster
localrules: pypsa1_retrieve_natura_raster, pypsa2_retrieve_cost_data, pypsa2_retrieve_natura_raster
localrules: pypsa1_download_copernicus_land_cover, pypsa1_build_powerplants, pypsa2_build_powerplants


rule build_load_data:
    message: "Use externally defined load time series in all scenarios."
    input: "data/electricity-demand-fully-electrified.csv"
    output: "pypsa-eur/resources/{scenario}/load.csv"
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
    output:
        s1 = f"build/results/scenarios/{SCENARIO_NAME_1}.nc",
        s2 = f"build/results/scenarios/{SCENARIO_NAME_2}.nc"
    run:
        copyfile(input.s1, output.s1)
        copyfile(input.s2, output.s2)


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
                spatial_res=config["resolution"]["space"],
                temporal_res=config["resolution"]["time"],
                opts="-" + opts),
            RESULT_PATH.format(
                scenario=SCENARIO_NAME_2,
                spatial_res=config["resolution"]["space"],
                temporal_res=config["resolution"]["time"],
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
    script: "../scripts/gsa/vis.py"
