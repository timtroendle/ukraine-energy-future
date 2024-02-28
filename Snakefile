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
        "build/results/gsa/capacities-power.csv",
        "build/results/gsa/capacities-energy.csv",
        "build/test-report.html",
        "build/test.success"


def pandoc_options(wildcards):
    suffix = wildcards["suffix"]
    if suffix == "html":
        return "--embed-resources --standalone --template template.html --to html5"
    elif suffix == "pdf":
        return "--template template.html --pdf-engine weasyprint"
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
        "report/cell.csl",
        "report/template.html",
        "report/more.css",
        "build/assumptions.csv",
        "build/gsa-parameters.csv",
        "data/final-energy.csv",
        "build/results/load.png",
        "build/results/load-validation.png",
        "build/results/lcoe-main.png",
        "build/results/capacities-main.png",
        "build/results/capacities-energy.csv",
        "build/results/generation-main.png",
        "build/results/time-series-main.png",
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
        "report/cell.csl",
        "report/supplementary.css",
        "report/template.html",
        "report/more.css",
        "build/results/lcoe-all.png",
        "build/results/capacities-all.png",
        "build/results/generation-all.png",
        "build/results/time-series-summer-zoom.png",
        "build/results/time-series-winter-zoom.png",
        "data/transmission-grid.png",
        "data/availability-offwind-ac.png",
        "data/availability-offwind-dc.png",
        "data/availability-onwind.png",
        "data/availability-solar.png",
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


rule push:
    message: "Package, zip, and move entire build."
    params:
        push_from_directory = config["push"]["from"],
        logs_directory = config["push"]["logs"],
        push_to_directory = config["push"]["to"]
    run:
        from datetime import datetime
        from pathlib import Path
        import shutil

        today = datetime.today().strftime('%Y-%m-%d')
        from_folder = Path(params.push_from_directory)
        logs_folder = Path(params.logs_directory)
        to_folder = Path(params.push_to_directory).expanduser()
        build_archive_filename = to_folder / "results" / f"ukraine-energy-future-{today}"
        logs_archive_filename = to_folder / "logs" / f"ukraine-energy-future-{today}"
        report_origin_pdf = from_folder / "report.pdf"
        report_destination_pdf = to_folder / "working-documents" / f"ukraine-energy-future-{today}.pdf"
        report_origin_docx = from_folder / "report.docx"
        report_destination_docx = to_folder / "working-documents" / f"ukraine-energy-future-{today}.docx"
        supplementary_origin = from_folder / "supplementary.pdf"
        supplementary_destination = to_folder / "working-documents" / f"ukraine-energy-future-{today}-supplementary.pdf"

        shutil.make_archive(build_archive_filename, 'zip', from_folder)
        shutil.make_archive(logs_archive_filename, 'zip', logs_folder)
        shutil.copyfile(report_origin_pdf, report_destination_pdf)
        shutil.copyfile(report_origin_docx, report_destination_docx)
        shutil.copyfile(supplementary_origin, supplementary_destination)


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
        main_scenarios = rules.run_scenarios.output.scenarios,
        main_scenario_logs = rules.run_scenarios.output.logs,
        gsa_scenarios = gsa_runs,
        gsa_scenario_logs = gsa_run_logs
    params:
        biomass_potential = config["pypsa-eur"]["biomass"]["potential"]
    resources:
        runtime = 60,
        mem_mb = 16000
    log: "build/test-report.html"
    output: "build/test.success"
    conda: "./envs/test.yaml"
    script: "./tests/test_runner.py"
