import pandas as pd
import altair as alt


DARK_GREY = "#424242"
TWO_COLUMN_WIDTH_PER_DISPLAY = 163
SINGLE_COLUMN_WIDTH = 473
DEFAULT_STEP = 14
SMALL_STEP = 10


def plot_all_capacities(power_capacities: pd.DataFrame, energy_capacities: pd.DataFrame, pre_war_capacities: dict,
                        nice_tech_names: dict[str: str], nice_scenario_names: dict[str, str],
                        scenarios: list[str], scenario_colors: dict[str: str],
                        generation_techs: list[str], energy_techs: list[str],
                        storage_techs: list[str], demand_techs: list[str]) -> alt.Chart:
    power_capacities = preprocess_capacities(power_capacities, pre_war_capacities, nice_tech_names, nice_scenario_names)
    energy_capacities = preprocess_capacities(energy_capacities, pre_war_capacities, nice_tech_names, nice_scenario_names).div(1000)
    capacities = pd.concat([power_capacities, energy_capacities])
    scenarios = list(map(nice_scenario_names.get, scenarios))
    scenario_colors = {nice_scenario_names.get(s, s): c for s, c in scenario_colors.items()}
    generation_techs = list(map(nice_tech_names.get, generation_techs))
    storage_techs = list(map(nice_tech_names.get, storage_techs))
    demand_techs = list(map(nice_tech_names.get, demand_techs))
    energy_techs = list(map(nice_tech_names.get, energy_techs))
    step = DEFAULT_STEP if len(scenarios) <= 3 else SMALL_STEP

    base = (
        alt
        .Chart(capacities.reset_index(), width=TWO_COLUMN_WIDTH_PER_DISPLAY)
    )

    panel_a = plot_capacities(base, generation_techs, "A", scenarios, scenario_colors, step=step)
    panel_b = plot_capacities(base, storage_techs, "B", scenarios, scenario_colors, step=step)
    panel_c = plot_capacities(
        base, energy_techs, "C", scenarios, scenario_colors, step,
        capacity_axis_label="Capacity (TWh)"
    )
    panel_d = plot_capacities(
        base, demand_techs, "D", scenarios, scenario_colors, step,
        carrier_axis_label=None, capacity_axis_label="Demand (GW)"
    )
    chart = (panel_a & panel_d) | (panel_b & panel_c)
    return configure_chart(chart)


def configure_chart(chart: alt.Chart) -> alt.Chart:
    return (
        chart
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom", columns=3, labelLimit=180)
    )


def plot_capacities(base: alt.Chart, techs: list[str], title: str, scenarios: list[str],
                    scenario_colors: dict[str: str], step: int, capacity_axis_label="Capacity (GW)",
                    carrier_axis_label="Technology") -> alt.Chart:
    return (
        base
        .transform_filter(alt.FieldOneOfPredicate(field='carrier', oneOf=techs))
        .transform_filter(alt.FieldOneOfPredicate(field='scenario', oneOf=scenarios))
        .encode(
            y=alt.Y("carrier:N").title(carrier_axis_label).sort(techs),
            x=alt.X("capacity:Q").title(capacity_axis_label).axis(labelFlush=False),
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
        .properties(title=title, height=alt.Step(step))
    )


def preprocess_capacities(sim_capacities: pd.DataFrame, pre_war_capacities: dict,
                          nice_tech_names: dict[str: str], nice_scenario_names: dict[str, str]) -> pd.DataFrame:
    pre_war = (
        pd.
        Series(pre_war_capacities)
        .rename("capacity")
        .rename_axis(index="carrier")
        .to_frame()
        .assign(scenario="pre-war")
        .reset_index()
        .set_index(["scenario", "carrier"])
        .rename(index=nice_scenario_names, level=0)
        .rename(index=nice_tech_names, level=1)
    )
    capacities = (
        sim_capacities
        .pipe(aggregate_offwind)
        .unstack()
        .rename("capacity")
        .rename_axis(index=["scenario", "carrier"])
        .rename(index=nice_scenario_names, level=0)
        .rename(index=nice_tech_names, level=1)
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
    chart = plot_all_capacities(
        power_capacities=pd.read_csv(snakemake.input.power_capacities, index_col=0),
        energy_capacities=pd.read_csv(snakemake.input.energy_capacities, index_col=0),
        pre_war_capacities=snakemake.params.pre_war,
        nice_tech_names=snakemake.params.nice_tech_names,
        nice_scenario_names=snakemake.params.nice_scenario_names,
        scenarios=snakemake.params.scenarios,
        scenario_colors=snakemake.params.scenario_colors,
        generation_techs=snakemake.params.generation_techs,
        storage_techs=snakemake.params.storage_techs,
        demand_techs=snakemake.params.demand_techs,
        energy_techs=snakemake.params.energy_techs
    )
    chart.save(snakemake.output[0])
