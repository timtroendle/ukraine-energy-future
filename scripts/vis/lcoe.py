import pandas as pd
import altair as alt


DARK_GREY = "#424242"
WIDTH = 411


def plot_lcoes(lcoes: pd.DataFrame, component_map: dict[str, str], component_colors: dict[str: str],
               scenarios: list[str], nice_scenario_names: dict[str, str]) -> alt.Chart:
    scenarios = list(map(nice_scenario_names.get, scenarios))
    lcoes = (
        lcoes
        .groupby(by=component_map, axis=1)
        .sum()
        .unstack()
        .rename_axis(index=["technology", "scenario"])
        .rename(index=nice_scenario_names, level=1)
        .rename("LCOE")
        .reset_index()
    )
    return (
        alt
        .Chart(lcoes, width=WIDTH)
        .transform_filter(alt.FieldOneOfPredicate(field='scenario', oneOf=scenarios))
        .encode(
            x=alt.X("LCOE").title("LCOE (€/MWh)"),
            y=alt.Y("scenario").title("Scenario"),
            color=(
                alt
                .Color("technology")
                .title("Component")
                .sort(list(component_colors.keys()))
                .scale(domain=list(component_colors.keys()), range=list(component_colors.values()))
            ),
            order=alt.Order("color_technology_sort_index:Q").sort("ascending")
        )
        .mark_bar()
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom")
    )


if __name__ == "__main__":
    chart = plot_lcoes(
        lcoes=pd.read_csv(snakemake.input.lcoes, index_col=0),
        component_map=snakemake.params.component_map,
        component_colors=snakemake.params.component_colors,
        nice_scenario_names=snakemake.params.nice_scenario_names,
        scenarios=snakemake.params.scenarios
    )
    chart.save(snakemake.output[0])
