import altair as alt
import pandas as pd

DARK_GREY = "#424242"
WIDTH = 441


def plot_load_validation(df: pd.DataFrame) -> alt.Chart:
    df = df.assign(isCovid=lambda df: df.ENTSOE_Demand_COVID_Removed_GW.isna())

    base = (
        alt
        .Chart(df.reset_index(), width=WIDTH)
        .encode(
            x=alt.X("ENTSOE_Demand_GW").scale(zero=False).title("Measured daily demand (GW)").axis(labelFlush=False),
            y=alt.Y("Demand_ninja_GW").scale(zero=False).title("Estimated daily demand (GW)"),

        )
    )
    points = (
        base
        .encode(color=alt.Color("isCovid:N").title("COVID19"))
        .mark_point(size=20, filled=True)
    )

    line = (
        base
        .encode(x=alt.X("ENTSOE_Demand_GW"), y=alt.Y("ENTSOE_Demand_GW"))
        .mark_line(color="black", strokeDash=[4, 4])
    )

    return (
        (points + line)
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
    )


if __name__ == "__main__":
    chart = plot_load_validation(
        df=pd.read_csv(snakemake.input.data, index_col=0, parse_dates=True)
    )
    chart.save(snakemake.output[0])
