
import pandas as pd
import altair as alt

DARK_GREY = "#424242"
WIDTH = 500


def plot_load(load: pd.DataFrame) -> alt.Chart:
    daily = load.div(1e3).resample("D").mean()
    return (
        alt
        .Chart(daily.stack().rename_axis(["utc_timestamp", "scenario"]).rename("load").reset_index(), width=WIDTH)
        .encode(
            x=alt.X("utc_timestamp:T", title="Time").axis(format='%B'),
            y=alt.Y("load", title="Mean daily electricity demand (GW)"),
            color=alt.Color("scenario", title="Demand scenario")
        )
        .mark_line()
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
    )


if __name__ == "__main__":
    chart = plot_load(
        load=pd.DataFrame({
            "low": pd.read_csv(snakemake.input.load_low, index_col=0, parse_dates=True)["UA"],
            "high": pd.read_csv(snakemake.input.load_high, index_col=0, parse_dates=True)["UA"]
        })
    )
    chart.save(snakemake.output[0])
