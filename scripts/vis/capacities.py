import pandas as pd
import altair as alt

SCENARIO_ORDER = ["pre-war", "nuclear-and-renewables", "only-renewables"]
SCENARIO_COLORS = ["#424242", "#FABC3C", "#679436"]
GENERATION_TECH_ORDER = ["fossil", "hydro", "nuclear", "biomass", "onwind", "solar"]
STORAGE_TECH_ORDER = ["PHS", "H2 electrolysis", "H2 fuel cell", "battery charger", "battery discharger"]
DARK_GREY = "#424242"
WIDTH_PER_DISPLAY = 185


def plot_capacities(capacities: pd.DataFrame, pre_war_capacities: dict,
                    nice_tech_names: dict[str: str]) -> alt.Chart:
    capacities = preprocess_capacities(capacities, pre_war_capacities, nice_tech_names)

    base = (
        alt
        .Chart(capacities.reset_index(), width=WIDTH_PER_DISPLAY)
    )

    nice_generation_tech_order = [nice_tech_names.get(t, t) for t in GENERATION_TECH_ORDER]
    generation = (
        base
        .transform_filter(alt.FieldOneOfPredicate(field='carrier', oneOf=nice_generation_tech_order))
        .encode(
            y=alt.Y("carrier:N").title("Carrier").sort(nice_generation_tech_order),
            x=alt.X("capacity:Q").title("Capacity (GW)"),
            color=(
                alt
                .Color("scenario:N")
                .title("Scenario")
                .sort(SCENARIO_ORDER)
                .scale(domain=SCENARIO_ORDER, range=SCENARIO_COLORS)
            ),
            yOffset=alt.YOffset("scenario:N").sort(SCENARIO_ORDER)
        )
        .mark_bar()
        .properties(title="A")
    )

    nice_storage_tech_order = [nice_tech_names.get(t, t) for t in STORAGE_TECH_ORDER]
    storage = (
        base
        .transform_filter(alt.FieldOneOfPredicate(field='carrier', oneOf=nice_storage_tech_order))
        .encode(
            y=alt.Y("carrier:N").title(None).sort(nice_storage_tech_order),
            x=alt.X("capacity:Q").title("Capacity (GW)"),
            color=(
                alt
                .Color("scenario:N")
                .title("Scenario")
                .sort(SCENARIO_ORDER)
                .scale(domain=SCENARIO_ORDER, range=SCENARIO_COLORS)
            ),
            yOffset=alt.YOffset("scenario:N").sort(SCENARIO_ORDER)
        )
        .mark_bar()
        .properties(title="B")
    )

    return (
        (generation | storage)
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom")
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
        .unstack()
        .rename("capacity")
        .rename_axis(index=["scenario", "carrier"])
        .rename(index=nice_tech_names)
    )
    return pd.concat([capacities.to_frame(), pre_war])


if __name__ == "__main__":
    chart = plot_capacities(
        capacities=pd.read_csv(snakemake.input.capacities, index_col=0),
        pre_war_capacities=snakemake.params.pre_war,
        nice_tech_names=snakemake.params.nice_tech_names
    )
    chart.save(snakemake.output[0])
