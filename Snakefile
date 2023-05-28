from snakemake.utils import min_version
PANDOC = "pandoc --filter pantable --filter pandoc-crossref --citeproc"

configfile: "config/default.yaml"
include: "./rules/sync.smk"
include: "./rules/pypsa.smk"
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
        "build/report.html",
        "build/test-report.html",
        "build/results/nuclear-and-renewables/elec_s_6_ec_lvopt_24H.nc",
        "build/results/only-renewables/elec_s_6_ec_lvopt_24H.nc"


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
        {PANDOC} report.md  --metadata-file=pandoc-metadata.yaml {params.options} \
        -o ../build/report.{wildcards.suffix}
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
         print("PyPSA-Eur data has not been cleaned.") # FIXME clean PyPSA data


rule test:
    conda: "envs/test.yaml"
    output: "build/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
