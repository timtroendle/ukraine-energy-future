all_scenarios = expand("build/results/scenarios/{scenario}.nc", scenario=config["scenarios"])


rule lcoe:
    message: "Calculate LCOE of all scenarios."
    input:
        scenarios = all_scenarios
    output:
        "build/results/lcoe.csv"
    params:
        ignore_existing = False
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


rule assumptions:
    message: "Generate table of assumptions."
    input:
        cost = "pypsa-eur/resources/only-renewables/costs.csv"
    params:
        technologies = [
            "onwind", "solar", "coal", "nuclear", "hydro", "PHS", "biomass",
            "battery inverter", "battery storage", "electrolysis", "fuel cell"
        ],
        parameters = ["FOM", "VOM", "investment", "lifetime"],
        parameter_to_format_string = {"FOM": ".1f", "investment": ".0f", "lifetime": ".0f"},
        biomass_parameters = config["pypsa-eur"]["biomass"]
    output:
        "build/assumptions.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/assumptions.py"


rule plot_load:
    message: "Plot load time series."
    input:
        load = "data/electricity-demand-fully-electrified.csv"
    output:
        "build/results/load.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/load.py"
