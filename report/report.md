==We are intending to publish this as a [*Report*](https://www.cell.com/joule/article-types) in Joule.==

# Highlights (expected)

* The economic potential for renewable energy is vastly larger than future demand.
* System cost of fully-renewable energy system close to a system with current share of nuclear.
* Autarkic system requires substantial storage, likely much less when connected to EUR (which we do not model).

# Introduction

Russia's invasion of Ukraine has led to massive damages to the national energy supply. Among the entire energy system, electricity supply has suffered the most: About two thirds of pre-war generation capacity is today either damaged, completely destroyed, or temporarily occupied [@UNDP:2023]. Even if sufficient generation capacity was available, it could not supply Ukraine's population and industry as almost half of transmission grid substations in government-controlled areas are out of operation following missile and drone attacks [@UNDP:2023].

Post-war recovery will require substantial investments into rebuilding Ukraine's energy supply. Apart from restoring energy access for the population as fast as possible, recovery should take energy security and climate objectives into consideration. Both require a reduction of the dependence of fossil fuels, a shift to low-carbon technologies, and large-scale electrification of the heat, mobility, and industrial sectors. Recovery of the electricity supply will therefore include not only rebuilding generation capacities to pre-war levels but an expansion to meet increased future demands.

Before the 2022 invasion, Ukraine's energy system was characterised by its natural resources and geopolitical situation. About 30% of its energy came from natural gas which was mainly imported despite the country's vast domestic reserves. Another roughly equally sized share came from low-carbon sources with nuclear being the largest contributor. Renewable energy played a minor role in Ukraine's energy mix but its importance was increasing. In particular solar photovoltaics was thriving, with a growth of about 1\ GW per year on average. Supported by international partners such as the OECD, the EU, EBRD and the IMF, Ukraine was exploring ways to develop a greener energy sector, notably through solar and wind power [@OECD:2021a]. The reforms implemented in the energy sector [@VerhovnaRadaUkrainy:2017; @VerhovnaRadaUkrainy:2021] (==Olena, can you send me translations of these titles?==) aimed to improve the investment climate and attract additional foreign capital.

It its 2021 nationally determined contributions, Ukraine committed to reaching carbon-neutrality in 2060 and to decreasing greenhouse gas emissions by 65% by 2030 [@Ukraine:2021]. During war-time, the Ukraine government has repeated its intention to adhere to these commitments (==source?==) and provided more detailed goals and instruments to reach them. By 2030, renewable energy should make up 27% of electricity generated in Ukraine, requiring an expansion of wind generation capacity of 10\ GW, up from ~1\ GW pre-war. A set of policy instruments is in preparation to enable this, including contracts for difference, and the facilitation of self consumption [@Ukraine:2023]. Natural gas is foreseen to be replaced by biomethane and hydrogen [@Directorate-GeneralforEnergy:2023], and the generation of nuclear power is intended to increase with the construction of 20 small modular nuclear reactors with a total capacity of 3.2\ GW [@Energoatom:2023].

Renewable energy may prosper in the future: With strong winds and solar irradiation and with high availability of eligible land, generation potentials of wind and solar power are high [@Kudria:2021;@Sukurova:2023]. Through the 2021 synchronisation of the European continental transmission grid and initiatives to foster trade with the European Union [@Directorate-GeneralforEnergy:2023], Ukraine may have the potential to become an exporter of renewable energy.

Here, we assess options for Ukraine to decarbonise its energy system by 2060. We apply a systems approach using a model that is highly resolved in space and time and that includes all energy sectors.

# Results

System cost of renewable only scenario with fully-electrified heat and transport sectors and with average cost and medium-high economic growth assumption is low (<70\ €/MWh, @fig:lcoe). Keeping nuclear at pre-war's levels shows 18% (12.2\ €/MWh) higher cost (@fig:lcoe). System cost includes generation, transmission, and storage capacities.

![**Cost breakdown of levelised system cost of electricity.** Cost for scenarios with fully-electrified heat and transport sectors, average cost assumptions, and medium-high economic growth assumption. See Supplemental Figure\ S1 for results based on both economic growth assumptions. Storage includes hydrogen, battery, and pumped-hydro storage.](build/results/lcoe-main.png){#fig:lcoe}

Meeting future electricity demand requires a substantial expansion of pre-war generation capacities in both scenarios (@fig:capacities-power), but especially in the renewables only scenario. This is explained by a switch to technologies with lower capacity factors (wind and solar), growing final energy demand based on population and economic growth, and an electrification of the heat and transport sectors (electricity demand ~5 times higher than pre-war).

![**Installed capacities and electricity demand.** Installed **(A)** generation and **(B)** storage capacities, and **(C)** average electricity demand of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors and with average cost assumptions and medium-high economic growth assumption. See Supplemental Figure\ S2 for results based on both economic growth assumptions.](build/results/capacities-power-main.png){#fig:capacities-power}

A fully-autarkic energy system requires substantial amounts of storage capacities in both scenarios, but slightly more in the renewables only scenario (@tbl:capacities-energy).

```table
---
caption: 'Installed storage capacities (GWh). {#tbl:capacities-energy}'
alignment: LRR
include: build/results/capacities-energy.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.3
    - 0.35
    - 0.35
---
```

In the only renewables scenario, biomass plays an important role (@fig:generation). Like onshore wind and solar, biomass generates roughly a third of the country's fully-electrified demand. In fact, all available biomass (513\ GWh) is used.

![**Electricity mix.** Electricity mix of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors and with average cost assumptions and medium-high economic growth assumption. See Supplemental Figure\ S3 for both economic growth assumptions.](build/results/generation-main.png){#fig:generation}

The cost of nuclear has the greatest impact on the cost penalty of the nuclear scenario (@fig:sensitivities). When Nuclear capital cost is doubled, the nuclear-and-renewables scenario becomes on average more than 90\ €/MWh more expensive compared with the only-renewables scenario. Doubling the capital cost of solar, wind, and biomass increases the relative cost of the only-renewables scenario, but to a lower degree (20--30\ €/MWh on average). Capital costs of storage technologies and marginal costs have low or almost no impact.

![**Sensitivities of cost penalty of nuclear scenario to model parameters.** Thick bars show the average increase in the cost penalty of the nuclear-and-renewables scenario, given a doubling of parameter values. Thin bars show the standard deviation of these effects and are therefore a proxy for the interactivity of the parameters: the larger the standard deviation, the more does the impact of the parameter depend on other parameter values.](build/results/gsa/sensitivities-annotated.png){#fig:sensitivities}

# Discussion

Both scenarios require large amounts of new generation and storage capacities. Given that biomass fuel is currently not constrained in the model, we can assume that storage capacities will increase even further once biomass fuel is constrained. However, we can also assume that storage capacities can decrease strongly with connecting Ukraine to neighbouring countries.

While required generation and storage capacities are large, there is no reason to believe these capacities could not be built. Ukraine can be powered by only renewables in the future at cost similar to a energy supply including nuclear.

---

The existing energy system of Ukraine has deep problems, caused primarily by severe damage to infrastructure facilities as a result of targeted missile attacks on energy infrastructure and equipment wear and tear, and its restoration requires significant investment, which may be comparable to the investment costs of building a new energy system. According to the calculations, the capital and operating costs of nuclear generation are 3-4 times higher than those of renewable energy sources. However, nuclear power has advantages, as it does not depend on seasonal fluctuations. The calculated need for storage capacities under both scenarios does not differ significantly or affect the choice.

The projected transition to an energy model based on the elimination of fossil fuels will be able to fully meet the growing demand and full electrification of industrial and residential energy consumption, including mobility, according to various forecasts of economic growth in Ukraine. The scenario with a combination of nuclear and green generation is more expensive, as it involves higher capital costs and LCOE, respectively. This also finds it prove in model-based analysis of 19 small modular reactors conducted by Steigerwald et al. The authors state that the costs to generate electricity is non-competitive when compared to current costs for generating electricity from renewable energy sources [@Steigerwald:2023].

Both scenarios envisage not only the generation of the required amount of electricity but also the decentralisation of the energy system, which is an essential condition for ensuring the country's energy security, given the geopolitical situation and the direct threat from Russia, which will potentially persist even after the war ends or enters a frozen phase. Although decentralisation of the energy sector will require the creation and maintenance of significant storage capacities, according to the calculations, this will not significantly affect the relative cost of electricity under both scenarios and, in the long run, may lead to a reduction in overall costs.

Modelling of a scenario using only renewable energy sources revealed a high potential of bioenergy, which could be equal to the pre-war needs for imports of fossil fuels. Given that biomass fuel is currently not constrained in the model, we can assume that storage capacities will increase even further once biomass fuel is constrained. However, we can also assume that storage capacities can decrease strongly with connecting Ukraine to neighbouring countries.

The potential of renewable energy in Ukraine far exceeds the estimated demand by 2060, which indicates that Ukraine has the potential to it can sell surplus energy generated. By using opportunities tendered by the Ukraine-EU electricity interconnection project, "European Network of Transmission System Operators for Electricity (ENTSO-E) of Ukraine," as an emergency response, Ukraine can transit from being a net importer of energy resources to an energy exporter in the post-war time. Ukraine's geographic location places it at the crossroads of major energy transit routes, including pipelines and transportation networks. This strategic position allows Ukraine to serve as a transit country for energy resources moving between Europe and Asia. Exports of energy resources, in turn, will help reduce the potential cost of electricity for domestic consumers, reducing the need to build storage facilities.

The development of new, cleaner energy systems will help reduce GHG emissions and improve environmental quality. Moreover, synchronized efforts of Ukraine and Europe will contribute to the achievement of the climate goals of both Ukraine and the European Union.

# Experimental procedures

## Resource Availability

### Lead Contact

Requests for further information, resources, and materials should be directed to the lead contact, Tim Tröndle (tim.troendle@usys.ethz.ch).

### Materials Availability

The study did not generate new materials.

### Data and Code Availability

The datasets generated during this study are available on Zenodo with DOI ==(after peer-review)==.

The model code and all analysis steps are publicly available as a reproducible Snakemake [@Koster:2012] workflow on Zenodo with DOI ==(after peer-review)==.

## Model overview

We model the energy system of Ukraine by soft-linking a demand projection and a cost-minimising capacity expansion model. We run two different scenarios and analyse the sensitivity of the cost difference to input parameters using a global sensitivity analysis. We describe each element of our analysis in the following.

## Energy demand projection

We project Ukraine's energy demand for the year 2060 in which Ukraine intends to be climate-neutral. We consider full electrification of the heat and transport sectors by this time as one option for complete decarbonisation. We generate an hourly time series of electricity demand in 2060 through soft-linking DESSTINEE [@Oreggioni:2022] which projects annual energy demands per sector and per energy carrier, with the Demand.ninja model [@Staffell:2023a] which upscales these demands to hourly resolution based on the prevailing weather. We use the resulting profile as an exogenous input to the capacity expansion model, ignoring demand flexibility.

Annual energy demands are modelled using macroeconomic projections informing changes in energy service demands. These changes are converted to final energy consumption using assumptions on fuel switching (towards electrification) and device efficiency. We assume a 1.0% annual reduction in population between 2020--2060  [@UnitedNationsDepartmentofEconomicandSocialAffairsPopulationDivision:2022] leading to a reduction from 43.9 million people in 2020 to 29.8 million in 2060. We assume a 4.3% annual growth in GDP per capita over the same period, taken from the high/high scenario of [@Smits:2019]. This sees GDP per capita growing from $13,000 in 2020 to $70,100 in 2060 (in Purchasing Power Parity terms). As the World Bank present a range of scenarios, we also test their low/high scenario, with 1.0% annual growth that gives $19,100 GDP per capita in 2060. Under these scenarios, Ukraine’s national GDP at PPP changes from $570\ bn in 2020 to $2,100\ bn in 2060 (high scenario) or remains flat at $570\ bn in 2060 (low scenario) as higher incomes are countered by population degrowth. We take the GDP split across economic sectors in 2020 from [@WorldBankGroup:2023; @WorldBankGroup:2023a]: 70% services, 20% industry, 10% agriculture. This sectoral split is modelled to follow the trend observed across European countries, in which the service sectors grow faster than other sectors as countries get richer [@Bossmann:2015]. This yields a GDP split in 2060 of 83% services, 13% industry and 4% agriculture in the high-growth scenario, and 74% services, 18% industry and 8% agriculture in the low-growth scenario.

DESSTINEE was updated to use 2019 as the baseline calibration year, as a proxy for 2020 which we avoid due to COVID impacts on energy demand.  Final energy demands were taken from the International Energy Agency [@IEA:2023] and ==(????????)==, shown in @tbl:final-energy-demand.  The split between end uses (for heating and cooling) is an output of the DESSTINEE model, derived from annual heating and cooling degree days.

The 2060 demand projection was constructed using a bottom-up approach in DESSTINEE.  The numerous assumptions for consumer behaviour, technology mix and efficiency for each sector were derived from the default values in DESSTINEE for Ukraine’s neighbouring European countries (Poland, Slovakia, Hungary and Romania), except for assumptions on the fuel basket for each sector which were overridden to give complete electrification (in line with the study’s aims).  Influential assumptions included the average efficiency of heat pumps (Coefficient of Performance increasing from 3 in 2020 to 4.25 in 2060) and the specific energy consumption of private cars decreasing from 1.1 to 0.3\ MJ/passenger-km from 2020 to 2060 (due to efficiency gains from electric vehicles).  While these assumptions are particularly uncertain (as are any technical projections for several decades into the future), we find the results are relatively insensitive to the overall scale of demand, and thus to these particular assumptions.

In @tbl:final-energy-demand, electricity consumption is seen to increase in all sectors in both the high and low-growth scenarios.  Total energy consumption also rises in the high-growth scenario, but falls in the low-growth scenario as increasing efficiency outweighs the increase in service demands due to greater prosperity.  The residential sector is the exception to this, due to the predominance of heating demand which is modelled to see very large efficiency gains due to the uptake of heat pumps.  Total demand for space heating sees a large reduction between 2020 and 2060 due to both this increasing efficiency and also the modelled increase in temperatures due to climate change: annual heating degree days are projected to fall by 17% from 3,150 to 2600 [@Staffell:2023a].  Conversely, demand for space cooling is projected to rise sharply for two reasons.  Firstly, the existing service is already provided by air conditioners, so efficiency improvements will be marginal technology improvements, rather than fuel switching from coal and oil. Secondly, climate change increases summer temperatures and thus service demands, with annual cooling degree days growing 35% from 580 to 780 [@Staffell:2023a].

```table
---
caption: 'Final energy demand by sector and carrier (TWh/yr). {#tbl:final-energy-demand}'
alignment: LRRRR
include: ../data/final-energy.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.2
    - 0.2
    - 0.2
    - 0.2
---
```

To derive a hourly load time series of the fully-electrified energy demand, in 2060, we disaggregate each sector’s annual energy demand over time. We distinguish between weather-invariant demands (industry, transport, and some end uses in residential and commercial buildings) and weather-dependent demands.

Weather-invariant demands are modelled within DESSTINEE using the partial decomposition approach outlined in [@Bossmann:2015].  Diurnal profiles (covering 24 hours) are defined for each sector and end use, for weekdays and weekends, summers and winters, and applied to the demand.  For example, heavy-industry is modelled as a constant 24/7 load, residential appliances follow a general profile that mirrors household occupation.  Electric vehicles are modelled as an equal mix of time-of-use charging (daytime and evenings) and smart charging (overnight).

The Demand.ninja model was used to model weather-dependent demands for space heating and cooling in residential and commercial buildings [@Staffell:2023a].  This derives a composite temperature index from four variables which influence thermal comfort: temperature, solar irradiance, wind speeds and humidity.  From this index, daily heating and cooling degree days were derived, which are then used to apportion heating and cooling demand over time. Open-access profiles were used, derived from the model’s global average parameters, including heating and cooling thresholds of 14*C and 20*C respectively [@Demand.ninja:2023]. Demands were then upscaled from daily to hourly resolution using the diurnal profiles for heating and cooling that were derived for Ukraine’s neighbouring countries [@Staffell:2023a].

This process was first validated against historical national electricity demand in Ukraine [@ENTSO-E:2023], by applying weather data from 2018 to 2020 to the model’s baseline year (@tbl:final-energy-demand). The first ten months of 2020 were removed from comparison due to the effects of the COVID-19 pandemic on human activity and thus on energy demand [@Mehlig:2021].  This validation (@fig:load-validation) revealed a strong correlation between modelled and actual demand, with an R2\ =\ 0.95 and residual standard error of +/-\ 0.5\ GW compared to mean demand of 17\ GW.

![**Estimated versus measured demand.** Demand is shown as daily averages for the years 2018--2021. In our validation, we exclude data from the COVID19 period January--October 2020 (highlighted) due to its irregular pattern.](build/results/load-validation.png){#fig:load-validation}

This same process was then applied to the 2060 scenarios, using the weather years of 2010--2014 to correspond with the years used for the wind and solar profiles (@fig:load). This profile was used in the capacity expansion model.

![**Electricity load projection used in the model.** Load is projected for the year 2060 and it is assumed that heat and transport are fully electrified. Horizontal axis shows the historical weather years used to project 2060 demand.](build/results/load.png){#fig:load}

## Capacity expansion model

To derive cost-minimal generation, storage, and transmission capacities, we apply PyPSA-Eur [@Horsch:2018] as capacity expansion model. Given the long-term perspective in our analysis, we apply a greenfield approach in which we do not consider any existing capacities other than the transmission grid and hydro generation capacities. We do not consider any connections to other countries, but model Ukraine as a stand-alone system.

The model has the option to expand generation and storage capacities of onshore wind, solar power, biomass, nuclear power, lithium-ion batteries, and hydrogen storage to meet demand in every hour of the year. Total system cost is derived by summing up annuities of investment, operation, and maintenance cost of all installed capacities (@tbl:technology-cost). The depreciation rate is set to 10% [@Andersson:2020]. The model finds the set of installed capacities with minimal total system cost.

The potential generation of solar and wind power is taken from PyPSA-Eur and has been derived in the following way.

First, the eligible area for wind and solar development is calculated for each region using *atlite* at 100\ m grid resolution [@Hofmann:2021].
For all renewable technologies, natural protection areas are excluded based on the World Database on Protected Areas (WDPA) [@UNEP-WCMC:2018].
Based on the Copernicus Global Land Cover dataset [@Buchhorn:2020], shrubland, herbaceous and sparse vegetation, and cropland are assumed to be eligible for wind and solar development.
In addition, while built-up areas are included for solar PV potentials, a distance of 1000\ m from built-up areas has to be kept for onshore wind turbines.
Offshore wind development is allowed up to a water depth of 50\ m, which is determined based on the GEBCO bathymetry dataset [@GEBCO:2015].
Furthermore, dense shipping lanes are excluded based on the World Bank's Global Shipping Traffic Density dataset [@Cerdeiro:2021].
Wind parks further out than 30\ km from shore are assumed to be DC-connected, whereas near-shore wind parks are assumed to be AC-connected.
For each renewable technology and region, the available area is multiplied with allowed deployment densities, approximating the socio-technical potential.
These densities are 3\ MW/km^2^ for onshore wind, 2\ MW/km^2^ for offshore wind, 1.7\ MW/km^2^ for solar.

Second, based on the eligible areas per region, capacity factor time series for wind and solar generation are calculated using *atlite*.
For this step, historical weather data for the years 2010--2014 from ECMWF's ERA5 reanalysis dataset [@Hersbach:2020] is used to convert wind speed and solar irradiance data to hourly capacity factors using models for typical wind turbines and solar panels.
The solar generation is calculated based on the incidence angle of solar irradiation, the panel tilt angle and the conversion efficiency of CdTe panel.
Power curves for a Vestas V112 3\ MW (onshore) and NREL 5\ MW (offshore) turbine are employed to map wind speeds scaled to hub height to power outputs.
The capacity factors of offshore wind generation are multiplied with a correction factor of 88.55% to approximately account for wake effects [@Bosch:2018].
Finally, the gridded dataset (0.25°x0.25°) is mapped onto the geographical shape of each region, using the available area as weighting.

We assume the future bioenergy potential to be 513\ TWh/yr (438\ TWh/yr solid biomass and 75\ TWh/yr biogas), following an assessment in [@Geletukha:2020]. Their projections include energy crops but exclude liquid fuels. Energy crops may require up to ~6.6% of total land.

```table
---
caption: 'Main technology assumptions of the model. Fixed operation and maintenance cost (FOM) is given in "% / year". Variable operation and maintenance cost (VOM) is given in "EUR / MWh". Investment cost are given in "EUR / kW". Technology lifetime is given in years. {#tbl:technology-cost}'
alignment: LRRRR
include: build/assumptions.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.2
    - 0.2
    - 0.2
    - 0.2
---
```

## Scenarios

We run two scenarios to assess the impact of nuclear power on system cost of a fully-decarbonised energy supply. In the first scenario, `only-renewables`, electricity is generated exclusively from renewable sources: solar, wind, hydro, and biomass. In the second scenario, `nuclear-and-renewables`, in addition to renewables, electricity can be generated from nuclear fuel. We set the minimal level of nuclear to its pre-war level of ~50% of electricity generation: ==60 GW (final number depends on energy demand)==.

## Global sensitivity analysis

To understand the impact of model parameters on cost difference between the two scenarios, we perform a global sensitivity analysis. We apply the Morris Method, an efficient method that estimates sensitivities accurately and efficiently [@Kristensen:2016; @Usher:2023]. We use the method to understand which model parameters impact the LCOE difference between the two scenarios the most. We assess eight cost parameters of the model (@tbl:gsa-parameters).

```table
---
caption: 'Parameter uncertainty ranges considered in the global sensitivity analysis. Minimum and maximum values are relative to the default values considered in the model. ==**UPDATE RANGES**== {#tbl:gsa-parameters}'
alignment: LRR
include: build/gsa-parameters.csv
include-encoding: UTF-8
markdown: True
---
```

The sampling of the Morris Method works in the following way: Given eight model parameters, the method chooses a random start value within the 8-dimensional space spanned by the model parameters. The method follows a one-step-at-a-time approach, in which each iteration of the method changes the value of a single parameter only. This approach leads to a trajectory of (8 + 1) points in the parameter space. To explore this space more comprehensively, we generate 15 such trajectories, creating 135 different combinations of parameter values. We run both scenarios for all combinations allowing us to attribute LCOE differences to parameter values. For each parameter, we observe its impact on LCOE difference once for each trajectory, a total of 15 times. We use the open-source Python library SALib to perform the global sensitivity analysis [@Herman:2017].

# References
