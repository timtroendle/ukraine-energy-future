---
subtitle: Supplementary material
abstract:
css:
    - https://timtroendle.github.io/signature-layout/layout/reset.css
    - https://timtroendle.github.io/signature-layout/layout/article.css
    - more.css
    - supplementary.css
figLabels: [S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15]
linkReferences: False
---

# @fig:lcoe: Levelised system cost of electricity for both economic growth scenarios

![**Cost breakdown of levelised system cost of electricity.** Cost for scenarios with fully-electrified heat and transport sectors, average cost assumptions, and low and high economic growth assumptions. Storage includes hydrogen, battery, and pumped-hydro storage.](build/results/lcoe-all.png){#fig:lcoe}

# @fig:generation: Electricity mix for both economic growth scenarios

![**Electricity mix.** Electricity mix of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors and with average cost assumptions and low and high economic growth assumption.](build/results/generation-all.png){#fig:generation}


# @fig:capacities: Installed power capacities for both economic growth scenarios

![**Installed capacities and electricity demand.** Installed **(A)** generation, **(B)** storage, and **(C)** energy storage capacities, and **(D)** average electricity demand of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors, average cost assumptions, and high and low economic growth assumption.](build/results/capacities-all.png){#fig:capacities}

# @fig:time-series-summer: Generation time series for a selected week in summer

![**Generation time series for the week starting 18 June 2012.** The areas show hourly generation of the different components of the energy system. In the low-biomass scenario, available biomass is constrained to 10% of the otherwise assumed value of 513\ TWh/yr. Storage includes hydrogen, battery, and pumped-hydro storage. Storage demand is not shown.](build/results/time-series-summer-zoom.png){#fig:time-series-summer}

# @fig:time-series-winter: Generation time series for a selected week in winter

![**Generation time series for the week starting 03 January 2011.** The areas show hourly generation of the different components of the energy system. In the low-biomass scenario, available biomass is constrained to 10% of the otherwise assumed value of 513\ TWh/yr. Storage includes hydrogen, battery, and pumped-hydro storage. Storage demand is not shown.](build/results/time-series-winter-zoom.png){#fig:time-series-winter}

# @fig:population: Population projection

![**Projection of Ukraine's population.** Source: UN World Population Prospects 2022 [@UnitedNationsDepartmentofEconomicandSocialAffairsPopulationDivision:2022]](build/results/population.png){#fig:population}

# @fig:transmission: Transmission grid

![**Transmission grid.** Pre-war transmission grid as used in the capacity expansion model.](../data/transmission-grid.png){#fig:transmission}

# @fig:solar: Eligible areas for solar deployment

![**Eligible areas for solar deployment.** Green areas on this map of Ukraine show the areas in which the model is able to deploy solar power.](../data/availability-solar.png){#fig:solar}

# @fig:onwind: Eligible areas for onshore wind deployment

![**Eligible areas for onshore wind deployment.** Green areas on this map of Ukraine show the areas in which the model is able to deploy onshore wind turbines.](../data/availability-onwind.png){#fig:onwind}

# @fig:offwind-dc: Eligible areas for DC-connected offshore wind deployment

![**Eligible areas for DC-connected offshore wind deployment.** Green areas on this map of Ukraine's Exclusive Economic Zone show the areas in which the model is able to deploy DC-connected offshore wind turbines.](../data/availability-offwind-dc.png){#fig:offwind-dc}

# @fig:offwind-ac: Eligible areas for AC-connected offshore wind deployment

![**Eligible areas for AC-connected offshore wind deployment.** Green areas on this map of Ukraine's Exclusive Economic Zone show the areas in which the model is able to deploy AC-connected offshore wind turbines.](../data/availability-offwind-ac.png){#fig:offwind-ac}


# Note S1: Capital cost of new nuclear power capacities

For the capital cost of new nuclear power plants we assume 7,940\ EUR~2015~/kW in the default scenarios and apply a range of 1,985-11,910 EUR~2015~/kW in the global sensitivity analysis.
The default value is derived from the average cost range given for total capital cost (6,900--12,200\ USD~2019~/kW) in the Lazard report v13 [@Lazard:2019].
This is converted to EUR~2015~/kW by applying an inflation rate of 2% and a currency conversion of 0.9\ EUR/USD, yielding 7,940\ EUR~2015~/kW.
<!-- USD2019 / 1.02^4 * 0.9 = EUR2015 -->
Newer cost estimates from the updated Lazard report v16 [@Lazard:2023] based on the new Plant Vogtle in Georgia (USA) are higher (8,475--13,925\ USD~2023~/kW), yielding 8,603\ EUR~2015~/kW according to the same conversion approach.
NREL's Annual Technology Baseline (ATB) database (based on EIA Annual Energy Outlook 2023 [@U.S.EnergyInformationAdministration:2023]) lists more optimistic estimates for 2023, stating a total cost of 6,768\ EUR~2015~/kW (8,811\ USD~2022~/kW).
In the ATB database, the lowest cost projections for 2030 and 2050 results in 5,937\ EUR~2015~/kW and 5,122\ EUR~2015~/kW, respectively.

# Bibliography
