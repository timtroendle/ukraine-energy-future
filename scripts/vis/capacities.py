import pandas as pd
import altair as alt


DARK_GREY = "#424242"
TWO_COLUMN_WIDTH_PER_DISPLAY = 185
SINGLE_COLUMN_WIDTH = TWO_COLUMN_WIDTH_PER_DISPLAY * 2.55


def plot_power_capacities(capacities: pd.DataFrame, pre_war_capacities: dict,
                          nice_tech_names: dict[str: str], scenarios: list[str],
                          scenario_colors: dict[str: str], generation_techs: list[str],
                          storage_techs: list[str], demand_techs: list[str]) -> alt.Chart:
    capacities = preprocess_capacities(capacities, pre_war_capacities, nice_tech_names)

    base = (
        alt
        .Chart(capacities.reset_index(), width=TWO_COLUMN_WIDTH_PER_DISPLAY)
    )

    panel_a = plot_capacities(base, generation_techs, "A", nice_tech_names, scenarios, scenario_colors)
    panel_b = plot_capacities(base, storage_techs, "B", nice_tech_names, scenarios, scenario_colors)
    panel_c = plot_capacities(
        base, demand_techs, "C", nice_tech_names, scenarios, scenario_colors,
        carrier_axis_label=None, capacity_axis_label="Demand (GW)"
    )
    chart = (panel_a | (panel_b & panel_c))
    return configure_chart(chart)


def plot_storage_capacities(capacities: pd.DataFrame, pre_war_capacities: dict,
                           nice_tech_names: dict[str: str], scenarios: list[str],
                           scenario_colors: dict[str: str], storage_techs: list[str]) -> alt.Chart:
    capacities = preprocess_capacities(capacities, pre_war_capacities, nice_tech_names).div(1000) # to TWh

    base = (
        alt
        .Chart(capacities.reset_index(), width=SINGLE_COLUMN_WIDTH)
    )

    chart = plot_capacities(
        base,
        storage_techs,
        "",
        nice_tech_names,
        scenarios,
        scenario_colors,
        capacity_axis_label="Capacity (TWh)")
    return configure_chart(chart)


def configure_chart(chart: alt.Chart) -> alt.Chart:
    return (
        chart
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom")
    )


def plot_capacities(base: alt.Chart, techs: list[str], title: str, nice_tech_names: dict[str: str],
                    scenarios: list[str], scenario_colors: dict[str: str], capacity_axis_label="Capacity (GW)",
                    carrier_axis_label="Carrier") -> alt.Chart:
    nice_generation_tech_order = [nice_tech_names.get(t, t) for t in techs]
    return (
        base
        .transform_filter(alt.FieldOneOfPredicate(field='carrier', oneOf=nice_generation_tech_order))
        .transform_filter(alt.FieldOneOfPredicate(field='scenario', oneOf=scenarios))
        .encode(
            y=alt.Y("carrier:N").title(carrier_axis_label).sort(nice_generation_tech_order),
            x=alt.X("capacity:Q").title(capacity_axis_label),
            color=(
                alt
                .Color("scenario:N")
                .title("Scenario")
                .sort(scenarios)
                .scale(domain=scenarios, range=[scenario_colors[s] for s in scenarios])
            ),
            yOffset=alt.YOffset("scenario:N").sort(scenarios)
        )
        .mark_bar()
        .properties(title=title)
    )


def preprocess_capacities(sim_capacities: pd.DataFrame, pre_war_capacities: dict,
                          nice_tech_names: dict[str: str]) -> pd.DataFrame:
    pre_war = (
        pd.
        Series(pre_war_capacities)
        .rename("capacity")
        .rename_axis(index="carrier")
        .to_frame()
        .assign(scenario="pre-war")
        .reset_index()
        .set_index(["scenario", "carrier"])
        .rename(index=nice_tech_names)
    )
    capacities = (
        sim_capacities
        .pipe(aggregate_offwind)
        .unstack()
        .rename("capacity")
        .rename_axis(index=["scenario", "carrier"])
        .rename(index=nice_tech_names)
    )
    return pd.concat([capacities.to_frame(), pre_war])


def aggregate_offwind(df: pd.DataFrame) -> pd.DataFrame:
    if "offwind-dc" in df.T.columns:
        return (
            df
            .T
            .assign(offwind=lambda df: df["offwind-dc"] + df["offwind-ac"])
            .T
        )
    else:
        return df


if __name__ == "__main__":
    match snakemake.params.type:
        case "power":
            chart = plot_power_capacities(
                capacities=pd.read_csv(snakemake.input.capacities, index_col=0),
                pre_war_capacities=snakemake.params.pre_war,
                nice_tech_names=snakemake.params.nice_tech_names,
                scenarios=snakemake.params.scenarios,
                scenario_colors=snakemake.params.scenario_colors,
                generation_techs=snakemake.params.generation_techs,
                storage_techs=snakemake.params.storage_techs,
                demand_techs=snakemake.params.demand_techs
            )
        case "energy":
            chart = plot_storage_capacities(
                capacities=pd.read_csv(snakemake.input.capacities, index_col=0),
                pre_war_capacities=snakemake.params.pre_war,
                nice_tech_names=snakemake.params.nice_tech_names,
                scenarios=snakemake.params.scenarios,
                scenario_colors=snakemake.params.scenario_colors,
                storage_techs=snakemake.params.storage_techs,
            )
        case _:
            raise ValueError("Unknown chart type.")
    chart.save(snakemake.output[0])
