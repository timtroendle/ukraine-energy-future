all_scenarios = expand("build/results/scenarios/{scenario}.nc", scenario=config["scenarios"])


rule lcoe:
    message: "Calculate LCOE of all scenarios."
    input:
        scenarios = all_scenarios
    output:
        "build/results/lcoe.csv"
    params:
        ignore_existing = False,
        opts_out = False
    conda: "../envs/default.yaml"
    script: "../scripts/lcoe.py"


rule capacities:
    message: "Calculate capacities of all scenarios."
    input:
        scenarios = all_scenarios
    output:
        power = "build/results/capacities-power.csv",
        energy = "build/results/capacities-energy.csv",
    conda: "../envs/default.yaml"
    script: "../scripts/capacities.py"


rule generation:
    message: "Calculate carrier generation in all scenarios."
    input:
        scenarios = all_scenarios
    output:
        energy = "build/results/generation.csv",
    conda: "../envs/default.yaml"
    script: "../scripts/generation.py"


rule assumptions:
    message: "Generate table of assumptions."
    input:
        cost = "pypsa-eur/resources/only-renewables-high/costs.csv"
    params:
        technologies = [
            "onwind", "offwind", "solar", "nuclear", "hydro", "PHS", "biomass",
            "battery inverter", "battery storage", "electrolysis", "fuel cell"
        ],
        parameters = ["FOM", "VOM", "investment", "lifetime"],
        parameter_to_format_string = {"FOM": ".1f", "investment": ".0f", "lifetime": ".0f"},
        biomass_parameters = config["pypsa-eur"]["biomass"]
    output:
        "build/assumptions.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/assumptions.py"


rule gsa_parameters:
    message: "Generate table of global sensitivity analysis parameters."
    params:
        parameters = config["global-sensitivity-analysis"]["parameters"]
    output:
        "build/gsa-parameters.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/gsa/parameters.py"


rule plot_load:
    message: "Plot load time series."
    input:
        load = "data/electricity-demand-fully-electrified.nc",
    output:
        "build/results/load.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/load.py"


rule plot_generation_capacities:
    message: "Plot generation capacities."
    input:
        capacities = rules.capacities.output.power
    params:
        pre_war = config["report"]["pre-war"]["capacities"],
        nice_tech_names = config["report"]["nice-names"]["technology"],
        scenarios = lambda wildcards, output: config["report"]["scenario-sets"][wildcards.scenario_set],
        scenario_colors = config["report"]["scenario-colors"]
    output:
        "build/results/capacities-power-{scenario_set}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/capacities.py"


rule plot_generation:
    message: "Plot generation shares."
    input:
        generation = rules.generation.output.energy
    params:
        pre_war = config["report"]["pre-war"]["electricity-mix"],
        nice_tech_names = config["report"]["nice-names"]["technology"],
        scenarios = lambda wildcards, output: config["report"]["scenario-sets"][wildcards.scenario_set]
    output:
        "build/results/generation-{scenario_set}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/generation.py"
