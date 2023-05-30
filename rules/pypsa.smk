import copy
from shutil import copyfile

from snakemake.io import load_configfile
from snakemake.utils import update_config

base_pypsa_config = load_configfile("pypsa-eur/config.default.yaml")

SCENARIO_NAME_1 = "nuclear-and-renewables"
SCENARIO_NAME_2 = "only-renewables"

scenario1_config = copy.deepcopy(base_pypsa_config)
scenario2_config = copy.deepcopy(base_pypsa_config)
update_config(scenario1_config, config["scenarios"][SCENARIO_NAME_1])
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


use rule * from pypsa_scenario1 as pypsa1_*
use rule * from pypsa_scenario2 exclude retrieve_databundle, retrieve_cutout, download_copernicus_land_cover, retrieve_load_data as pypsa2_*


####################
# The following moves results from pypsa-eur to build folder of this project
####################


rule run_scenarios:
    input:
        s1 = "pypsa-eur/results/networks/nuclear-and-renewables/elec_s_6_ec_lvopt_24H.nc",
        s2 = "pypsa-eur/results/networks/only-renewables/elec_s_6_ec_lvopt_24H.nc"
    output:
        s1 = "build/results/scenarios/nuclear-and-renewables.nc",
        s2 = "build/results/scenarios/only-renewables.nc"
    run:
        copyfile(input.s1, output.s1)
        copyfile(input.s2, output.s2)
