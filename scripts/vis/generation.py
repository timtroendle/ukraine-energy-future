import pandas as pd
import altair as alt


GENERATION_TECH_ORDER = ["nuclear", "fossil", "hydro", "biomass", "offwind", "onwind", "solar"]
TECH_COLORS = ["#A01914", "#424242", "#4F6DB8", "#679436", "#6851A3", "#9483C1", "#FABC3C"]
DARK_GREY = "#424242"
WIDTH = 315


def plot_generation(generation: pd.DataFrame, pre_war_generation: dict,
                    nice_tech_names: dict[str: str], scenarios: list[str],
                    nice_scenario_names: dict[str, str]) -> alt.Chart:
    generation = preprocess_generation(generation, pre_war_generation, nice_tech_names, nice_scenario_names)
    scenarios = list(map(nice_scenario_names.get, scenarios))

    base = (
        alt
        .Chart(generation.reset_index(), width=WIDTH)
    )

    nice_generation_tech_order = [nice_tech_names.get(t, t) for t in GENERATION_TECH_ORDER]
    return (
        base
        .transform_filter(alt.FieldOneOfPredicate(field='carrier', oneOf=nice_generation_tech_order))
        .transform_filter(alt.FieldOneOfPredicate(field='scenario', oneOf=scenarios))
        .encode(
            color=(
                alt
                .Color("carrier:N")
                .title("Technology")
                .scale(domain=nice_generation_tech_order, range=TECH_COLORS)
                .sort(nice_generation_tech_order)
            ),
            x=alt.X("generation:Q").title("Generation").stack("normalize").axis(labelFlush=False),
            y=alt.Y("scenario:N").title("Scenario").sort(scenarios),
            order=alt.Order("color_carrier_sort_index:Q").sort("ascending")
        )
        .mark_bar()
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY)
    )


def preprocess_generation(sim_generation: pd.DataFrame, pre_war_generation: dict,
                          nice_tech_names: dict[str: str], nice_scenario_names: dict[str, str]) -> pd.DataFrame:
    sim_generation = (
        sim_generation
        .T
        .assign(offwind=lambda df: df["offwind-dc"] + df["offwind-ac"])
        .T
        .unstack()
        .rename("generation")
        .rename_axis(index=["scenario", "carrier"])
        .rename(index=nice_scenario_names, level=0)
        .rename(index=nice_tech_names, level=1)
    )
    pre_war = (
        pd
        .Series(pre_war_generation)
        .rename("generation")
        .rename_axis(index="carrier")
        .rename(index=nice_tech_names)
        .to_frame()
        .assign(scenario="pre-war")
        .reset_index()
        .set_index(["scenario", "carrier"])
        .rename(index=nice_scenario_names, level=0)
    )
    return pd.concat([sim_generation.to_frame(), pre_war])


if __name__ == "__main__":
    chart = plot_generation(
        generation=pd.read_csv(snakemake.input.generation, index_col=0),
        pre_war_generation=snakemake.params.pre_war,
        nice_tech_names=snakemake.params.nice_tech_names,
        nice_scenario_names=snakemake.params.nice_scenario_names,
        scenarios=snakemake.params.scenarios
    )
    chart.save(snakemake.output[0])
