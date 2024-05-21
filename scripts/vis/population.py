import altair as alt
import pandas as pd

DARK_GREY = "#424242"
domain = [
    "Historical",
    "Low fertility scenario",
    "Medium fertility scenario",
    "High fertility scenario",
]
range_ = ["black", "#A0CB71", "#679436", "#A0CB71"]
WIDTH = 505


def plot_population(population: pd.DataFrame) -> alt.Chart:
    base = alt.Chart(population, width=WIDTH).encode(
        x=alt.X("Time:T").title(None),
        y=alt.Y("population").title("Population (millions)"),
        color=alt.Color("Variant").scale(domain=domain, range=range_).title(None),
    )

    medium = base.transform_filter(
        alt.FieldOneOfPredicate(
            field="Variant", oneOf=["Historical", "Medium fertility scenario"]
        )
    ).mark_line(clip=True, strokeWidth=3)
    high_low = base.transform_filter(
        alt.FieldOneOfPredicate(
            field="Variant", oneOf=["Low fertility scenario", "High fertility scenario"]
        )
    ).mark_line(clip=True)

    return (
        (high_low + medium)
        .configure(font="Lato")
        .configure_title(anchor="start", fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom")
    )


if __name__ == "__main__":
    chart = plot_population(population=pd.read_feather(snakemake.input[0]))
    chart.save(snakemake.output[0])
