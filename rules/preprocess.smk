localrules: download_population_data


rule download_population_data:
    message: "Get population data from UN."
    params:
        url = lambda wildcards: config["data-sources"][f"population-{wildcards.scenario}"]
    output: protected("data/automatic/population-{scenario}.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sSLo {output} '{params.url}'"


rule population_data:
    message: "Pre-process population data."
    input:
        medium = "data/automatic/population-medium.zip",
        high_low = "data/automatic/population-high-low.zip"
    output:
        "build/data/population.feather"
    script: "../scripts/preprocess/population.py"
