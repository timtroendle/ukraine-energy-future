
import pandas as pd
import altair as alt

DARK_GREY = "#424242"
WIDTH = 500


def plot_load(load: pd.Series) -> alt.Chart:
    daily = load.div(1e3).resample("D").mean()
    return (
        alt
        .Chart(daily.reset_index(), width=WIDTH)
        .encode(x=alt.X("utc_timestamp", title="Time"), y=alt.Y("UA", title="Mean daily electricity demand (GW)"))
        .mark_line()
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
    )


if __name__ == "__main__":
    chart = plot_load(
        load=pd.read_csv(snakemake.input.load, index_col=0, parse_dates=True)
    )
    chart.save(snakemake.output[0])
