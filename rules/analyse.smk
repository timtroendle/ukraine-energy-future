rule lcoe:
    message: "Calculate LCOE of all scenarios."
    input:
        scenarios = expand("build/results/scenarios/{scenario}.nc", scenario=config["scenarios"])
    output:
        "build/results/lcoe.csv"
    params:
        ignore_existing = False
    conda: "../envs/default.yaml"
    script: "../scripts/lcoe.py"

