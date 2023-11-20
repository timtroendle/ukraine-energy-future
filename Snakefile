from snakemake.utils import min_version
PANDOC = "pandoc --filter pantable --filter pandoc-crossref --citeproc"

configfile: "config/default.yaml"
include: "./rules/sync.smk"
include: "./rules/pypsa.smk"
include: "./rules/analyse.smk"
include: "./rules/utils.smk"
localrules: all, report, clean
min_version("7.8")

onstart:
    shell("mkdir -p build/logs")
onsuccess:
    if "email" in config.keys():
        shell("echo "" | mail -s 'ukraine-energy-future succeeded' {config[email]}")
onerror:
    if "email" in config.keys():
        shell("echo "" | mail -s 'ukraine-energy-future failed' {config[email]}")


rule all:
    message: "Run entire analysis and compile report."
    input:
        "build/report.docx",
        "build/report.pdf",
        "build/supplementary.pdf",
        "build/test-report.html",
        "build/test.success"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--embed-resources --standalone --to html5"
    elif suffix == "pdf":
        return "--pdf-engine weasyprint"
    elif suffix == "docx":
        return []
    else:
        raise ValueError(f"Cannot create report with suffix {suffix}.")


rule report:
    message: "Compile report.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/report.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "build/assumptions.csv",
        "build/gsa-parameters.csv",
        "data/final-energy.csv",
        "build/results/load.png",
        "build/results/load-validation.png",
        "build/results/lcoe.csv",
        "build/results/capacities-power-main.png",
        "build/results/capacities-energy.csv",
        "build/results/generation-main.png",
        "build/results/gsa/sensitivities-annotated.png",
    params: options = pandoc_options
    output: "build/report.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} report.md --metadata-file=pandoc-metadata.yaml {params.options} \
        -f markdown+mark \
        -o ../build/report.{wildcards.suffix}
        """


rule supplementary:
    message: "Compile supplementary.{wildcards.suffix}."
    input:
        "report/literature.yaml",
        "report/supplementary.md",
        "report/pandoc-metadata.yaml",
        "report/apa.csl",
        "report/supplementary.css",
        "build/results/capacities-power-all.png",
        "build/results/generation-all.png",
    params: options = pandoc_options
    output: "build/supplementary.{suffix}"
    wildcard_constraints:
        suffix = "((html)|(pdf)|(docx))"
    conda: "envs/report.yaml"
    shadow: "minimal"
    shell:
        """
        cd report
        ln -s ../build .
        {PANDOC} supplementary.md --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/supplementary.{wildcards.suffix}
        """


rule dag:
     message: "Plot dependency graph of the workflow."
     conda: "envs/dag.yaml"
     shell:
         """
         snakemake --rulegraph > build/dag.dot
         dot -Tpdf -o build/dag.pdf build/dag.dot
         """


rule clean: # removes all generated results
    message: "Remove all build results but keep downloaded data."
    run:
         import shutil

         shutil.rmtree("build")
         print("Data downloaded to data/ has not been cleaned.")
         print("PyPSA-Eur data has not been cleaned.")


rule test:
    message: "Run tests"
    input:
        test_dir = "tests",
        tests = map(str, Path("tests").glob("**/test_*.py")),
        scenarios = expand("build/results/scenarios/{scenario}.nc", scenario=config["scenarios"].keys())
    log: "build/test-report.html"
    output: "build/test.success"
    conda: "./envs/test.yaml"
    script: "./tests/test_runner.py"
