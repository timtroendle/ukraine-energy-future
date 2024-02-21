from dataclasses import dataclass, field, InitVar

import pandas as pd
import pypsa
import altair as alt

DARK_GREY = "#424242"
WIDTH = 151
COMPONENT_ORDER = {
    "Wind, water, sun": 0,
    "Storage": 1,
    "Biomass": 2,
    "Nuclear": 3
}

@dataclass
class Scenario:
    n: pypsa.Network
    name: str = field(init=False)
    nice_name: str = field(init=False)
    nice_names: InitVar[dict[str, str] | None] = None
    opt_in_name: InitVar[bool] = False

    def __post_init__(self, nice_names, opt_in_name):
        self.name = self.n._meta["run"]["name"]
        if nice_names is not None and self.name in nice_names.keys():
            self.nice_name = nice_names[self.name]
        else:
            self.nice_name = self.name
        if opt_in_name:
            self.name = self.name + "--" + "gsa"


def plot_time_series(scenarios: list[Scenario], start_date: str, end_date: str, resolution: str,
                     colors: dict[str: str]) -> alt.Chart:
    colors = {component: colors[component] for component in COMPONENT_ORDER.keys()}
    df = all_generation(scenarios, start_date, end_date, resolution)
    time_format = "%Y" if resolution=="M" else "%d %b %y"
    chart = (
        alt
        .Chart(df.reset_index(names="timestamp"), width=WIDTH, height=WIDTH)
        .encode(
            x=alt.X('timestamp').title(None).axis(format=time_format),
            y=alt.Y('generation').title("Generation (GW)"),
            color=(
                alt
                .Color('tech')
                .title("Component")
                .sort(list(COMPONENT_ORDER.keys()))
                .scale(domain=list(colors.keys()), range=list(colors.values()))
            ),
            order=alt.Order('order', sort='descending')
        )
        .mark_area(line=False)
        .facet("Scenario")
        .configure(font="Lato")
        .configure_title(anchor='start', fontSize=12, color=DARK_GREY)
        .configure_axis(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_header(titleColor=DARK_GREY, labelColor=DARK_GREY)
        .configure_legend(titleColor=DARK_GREY, labelColor=DARK_GREY, orient="bottom")
    )
    return chart


def all_generation(scenarios: list[Scenario], start_date: str, end_date: str, resolution: str) -> pd.DataFrame:
    return pd.concat([
        generation(s.n, start_date, end_date, resolution).assign(Scenario=s.nice_name)
        for s in scenarios
    ])


def generation(n: pypsa.Network, start_date: str, end_date: str, resolution: str) -> pd.DataFrame:
    gen = (
        pd
        .DataFrame({
            tech_name: n.generators_t.p[[g for g in n.generators.index if tech_name in g]].sum(axis=1)
            for tech_name in ["solar", "nuclear", "wind", "hydro"]
        })
    )
    storage = (
        pd
        .DataFrame({
            tech_name: n.links_t.p0[[g for g in n.links.index if tech_name in g]].sum(axis=1)
            for tech_name in ["biomass", "Fuel Cell", "battery discharger"]
        })
    )
    return (
        pd
        .concat([gen, storage], axis=1)
        .div(1000) # to GW
        .assign(
            storage=lambda df: df[["Fuel Cell", "battery discharger"]].sum(axis="columns"),
            wws=lambda df: df[["wind", "solar", "hydro"]].sum(axis="columns")
        )
        .drop(columns=["Fuel Cell", "battery discharger", "wind", "solar", "hydro"])
        .rename(columns={"wws": "Wind, water, sun"})
        .rename(columns=lambda name: name.capitalize())
        .loc[start_date:end_date]
        .resample(resolution)
        .mean()
        .melt(var_name="tech", value_name="generation", ignore_index=False)
        .assign(order=lambda df: df.tech.map(COMPONENT_ORDER))
    )

# TODO include PHS

if __name__ == "__main__":
    scenarios = (
        Scenario(n=pypsa.Network(path_to_n), nice_names=snakemake.params.nice_scenario_names, opt_in_name=False, )
        for path_to_n in snakemake.input.scenarios
    )
    scenarios = (s for s in scenarios if s.name in snakemake.params.scenarios)
    chart = plot_time_series(
        scenarios=scenarios,
        start_date=snakemake.params.start_date,
        end_date=snakemake.params.end_date,
        resolution=snakemake.params.resolution,
        colors=snakemake.params.colors
    )
    chart.save(snakemake.output[0])
