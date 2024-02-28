import pandas as pd
import altair as alt

DARK_GREY = "#424242"
WIDTH = 421


def visualise_sensitivities(sensitivities: pd.DataFrame, cost_penalty: float) -> alt.Chart:
    sensitivities = (
        sensitivities
        .assign(
            upper=lambda df: df.mu + df.sigma,
            lower=lambda df: df.mu - df.sigma,
            baseline=-cost_penalty,
            names=lambda df: (
                df
                .names
                .str
                .replace("onwind", "wind")
                .str
                .replace("biomassinflow", "Biomass")
                .str
                .replace("\+c", " capital cost", regex=True)
                .str
                .replace("\+m", " marginal cost", regex=True)
                .str
                .replace("\+l", " potential", regex=True)
                .str
                .capitalize()
            )
        )
    )
    chart = (
        alt
        .Chart(sensitivities, width=WIDTH)
        .encode(
            y=alt.Y("names", title="Parameter").sort(alt.EncodingSortField("mu_star", op="max", order="descending")),
            x=alt.X("mu").title("Impact on cost penalty of nuclear power (EUR/MWh)").axis(labelFlush=False)
        )
    )
    bar = chart.mark_bar()
    errorbar = (
        chart
        .mark_errorbar()
        .encode(
            x=alt.X("lower", title=""),
            x2=alt.X2("upper", title=None)
        )
    )
    baseline = (
        alt
        .Chart(sensitivities)
        .mark_rule(color=DARK_GREY, strokeDash=[2], opacity=0.8)
        .encode(x=alt.X("baseline:Q"))
    )

    return (
        (baseline + bar + errorbar)
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
    )


if __name__ == "__main__":
    lcoe = pd.read_csv(snakemake.input.lcoe, index_col=0)["System LCOE (EUR/MWh)"]
    cost_penalty = lcoe["nuclear-and-renewables-high"] - lcoe["only-renewables-high"]
    chart = visualise_sensitivities(
        pd.read_csv(snakemake.input.sensitivities),
        cost_penalty=cost_penalty
    )
    chart.save(snakemake.output[0])
