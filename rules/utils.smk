localrules: render_vega_lite_to_png

rule render_vega_lite_to_png:
    message: "Render Vega Lite spec {wildcards.filename}.json to png."
    input:
        json = "build/{path}/{filename}.vega.json"
    output:
        png = "build/{path}/{filename}.png"
    conda: "../envs/vega.yaml"
    # vl2png not usable because of https://github.com/queryverse/VegaLite.jl/issues/383
    shell: "vl2vg {input.json} | vg2png --scale 4 > {output.png}" # scale to >300dpi
