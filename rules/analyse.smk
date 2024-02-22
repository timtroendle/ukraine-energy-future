all_scenarios = expand("build/results/scenarios/{scenario}.nc", scenario=config["scenarios"])


rule lcoe:
    message: "Calculate LCOE of all scenarios."
    input:
        scenarios = all_scenarios
    output:
        "build/results/lcoe.csv"
    params:
        opts_out = False
    conda: "../envs/default.yaml"
    script: "../scripts/lcoe.py"


rule capacities:
    message: "Calculate capacities of all scenarios."
    input:
        scenarios = all_scenarios
    params:
        opts_out = False
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
        cost = "pypsa-eur/resources/only-renewables-high/costs.csv",
        sources = "data/sources-cost-data.csv"
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
    params:
        colors = {"2060 Low Growth": "#A0CB71", "2060 High Growth": "#679436"}
    output:
        "build/results/load.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/load.py"


rule plot_load_validation:
    message: "Plot validation of load time series."
    input:
        data = "data/ninja.historical.validation.csv",
    output:
        "build/results/load-validation.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/validation.py"


rule plot_power_capacities:
    message: "Plot power capacities."
    input:
        capacities = rules.capacities.output.power
    params:
        type = "power",
        pre_war = config["report"]["pre-war"]["capacities"],
        nice_tech_names = config["report"]["nice-names"]["technology"],
        nice_scenario_names = config["report"]["nice-names"]["scenario"],
        scenarios = lambda wildcards, output: config["report"]["scenario-sets"][wildcards.scenario_set],
        scenario_colors = config["report"]["scenario-colors"],
        generation_techs = ["fossil", "hydro", "nuclear", "biomass", "onwind", "offwind", "solar"],
        storage_techs = ["PHS", "H2 electrolysis", "H2 fuel cell", "battery charger", "battery discharger"],
        demand_techs = ["average-load"]
    output:
        "build/results/capacities-power-{scenario_set}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/capacities.py"


rule plot_energy_capacities:
    message: "Plot energy capacities."
    input:
        capacities = rules.capacities.output.energy
    params:
        type = "energy",
        pre_war = config["report"]["pre-war"]["capacities"],
        nice_tech_names = config["report"]["nice-names"]["technology"],
        nice_scenario_names = config["report"]["nice-names"]["scenario"],
        scenarios = lambda wildcards, output: config["report"]["scenario-sets"][wildcards.scenario_set],
        scenario_colors = config["report"]["scenario-colors"],
        storage_techs = ["battery", "H2"],
    output:
        "build/results/capacities-energy-{scenario_set}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/capacities.py"


rule plot_generation:
    message: "Plot generation shares."
    input:
        generation = rules.generation.output.energy
    params:
        pre_war = config["report"]["pre-war"]["electricity-mix"],
        nice_tech_names = config["report"]["nice-names"]["technology"],
        nice_scenario_names = config["report"]["nice-names"]["scenario"],
        scenarios = lambda wildcards, output: config["report"]["scenario-sets"][wildcards.scenario_set]
    output:
        "build/results/generation-{scenario_set}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/generation.py"


rule plot_lcoes:
    message: "Plot LCOE breakdown."
    input:
        lcoes = rules.lcoe.output[0]
    params:
        component_map = config["report"]["components"]["mapping"],
        component_colors = config["report"]["components"]["colors"],
        scenarios = lambda wildcards, output: config["report"]["scenario-sets"][wildcards.scenario_set],
        nice_scenario_names = config["report"]["nice-names"]["scenario"]
    output:
        "build/results/lcoe-{scenario_set}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/lcoe.py"


rule plot_time_series:
    message: "Plot generation time series."
    input:
        scenarios = all_scenarios
    params:
        scenarios = config["report"]["scenario-sets"]["time-series"],
        start_date = lambda wildcards, output: config["report"]["time-series-plots"][wildcards.plot]["start"],
        end_date = lambda wildcards, output: config["report"]["time-series-plots"][wildcards.plot]["end"],
        resolution = lambda wildcards, output: config["report"]["time-series-plots"][wildcards.plot]["resolution"],
        colors = config["report"]["components"]["colors"],
        nice_scenario_names = config["report"]["nice-names"]["scenario"]
    output:
        "build/results/time-series-{plot}.vega.json"
    conda: "../envs/default.yaml"
    script: "../scripts/vis/time_series.py"
