# Highlights

* The economic potential for renewable energy is vastly larger than future demand.
* System cost of fully-renewable energy supply is low.
* Both nuclear power and biomass-based generation are not essential.
* We couple open models and data in an automated workflow for reproducible research.

# Introduction

Russia's 2022 invasion of Ukraine has led to massive damage to the national energy supply. Among the energy system, electricity supply has suffered the most: One year into the war, about two thirds of pre-war generation capacity was either damaged, completely destroyed, or temporarily occupied [@UNDP:2023]. Even if sufficient generation capacity was available, it could not supply Ukraine's population and industry, as almost half of transmission grid substations in government-controlled areas are out of operation following missile and drone attacks [@UNDP:2023].

Post-war recovery will require substantial investments into rebuilding Ukraine's energy supply. Apart from restoring energy access for the population as fast as possible, recovery should take longer-term energy security and climate objectives into consideration. Both require reducing the dependence on fossil fuels, via shifting to low-carbon technologies, and large-scale electrification of the heat, mobility, and industrial sectors. Restoring the nation’s electricity supply will therefore include not only rebuilding generation capacities to pre-war levels but expanding them to meet increased future demands.

Before the 2022 invasion, Ukraine's energy system was characterised by its natural resources and geopolitical situation. About 30% of its energy came from natural gas, of  which about a third was imported [@EnergyInstitute:2023]. A similar share of energy came from low-carbon sources with nuclear being the largest contributor [@EnergyInstitute:2023]. Renewable energy from wind, water, and solar sources played a minor but growing role in Ukraine's energy mix. Solar photovoltaics was growing particularly fast, with an annual growth rate of about 36% [@EnergyInstitute:2023]. Supported by international partners such as the OECD, the EU, EBRD, and the IMF, Ukraine was exploring ways to develop a greener energy sector, notably through solar and wind power [@OECD:2021a]. The reforms implemented in the energy sector [@VerhovnaRadaUkrainy:2017; @VerhovnaRadaUkrainy:2021] aimed to improve the investment climate and attract additional foreign capital.

In its 2021 Nationally Determined Contribution, Ukraine committed to reaching carbon-neutrality in 2060 and to decreasing greenhouse gas emissions by 65% by 2030 [@Ukraine:2021]. During war-time, Ukraine’s government has reemphasised its intention to adhere to these commitments [@CabinetofMinistersofUkraine:2022] and provided more detailed goals and instruments to reach them. By 2030, renewable energy should make up 27% of electricity generated in Ukraine, requiring an expansion of wind generation capacity to a total of 10\ GW, up from ~1\ GW pre-war. A set of policy instruments is in preparation to enable this, including contracts for difference, and support of self consumption [@Ukraine:2023]. Natural gas is foreseen to be replaced by biomethane and hydrogen [@Directorate-GeneralforEnergy:2023], and the generation of nuclear power is intended to increase through the construction of 20 small modular nuclear reactors with a total capacity of 3.2\ GW [@Energoatom:2023].

Renewable energy may prosper in the future: With high wind speeds and solar irradiation and substantial availability of eligible land, generation potentials of wind and solar power are high [@Kudria:2021;@Sukurova:2023]. Through the 2021 synchronisation with the European continental transmission grid and initiatives to foster trade with the European Union [@Directorate-GeneralforEnergy:2023], Ukraine has the potential to become an exporter of renewable energy in the future [@Sukurova:2023].

However, despite high scientific and non-scientific attention towards rebuilding Ukraine and its energy supply [@Worldbank:2023; @Becker:2022; @IKEM:2023; @GreenDealUkraine:2024], Ukraine’s options to decarbonise its energy system by 2060 and their implications are not properly understood. Previous studies have shown that 1.5 °C mitigation pathways are possible at small or no cost penalty: @Chepeliev:2023 identified only a small cost penalty compared with a fossil-fuel based reference scenario using a substantial expansion of biomass-based generation, largely ignoring the vast potentials of wind and solar power. In contrast, solar power builds the foundation of European scenarios that include Ukraine in a study by @Breyer:2023, in which levelised cost of energy decreases compared with today’s levels. A related study applying the same model [@Child:2017a] finds decreasing cost even when Ukraine is modelled in isolation which leads to higher costs than when trade with neighbours is possible. Another study by @Osadcha:2024 calls for the expansion of nuclear power in order to replace fossil fuels and decrease import dependencies. None of these studies compares the different technology options and assesses the large uncertainties in future energy demand stemming from developments in population, economic growth, and technology in the coming decades.

Here, we systematically assess options for Ukraine to decarbonise its energy system by 2060 using detailed open datasets and models. We apply a systems approach using models that are highly resolved in space and time and that include all energy sectors. We combine a bottom-up energy demand model with a cost-minimising capacity expansion model and assess uncertainties using a global sensitivity analysis. We find that a full decarbonisation of the energy sector is possible at low cost, and that neither nuclear nor biomass are essential elements in the generation portfolio.

# Methods

We model the energy system of Ukraine by soft-linking a demand projection and a cost-minimising capacity expansion model. We run five different scenarios and analyse the sensitivity to input parameters using a global sensitivity analysis. We combine all elements using the workflow tool Snakemake [@Koster:2012] and we perform our calculations on ETH Zurich’s high-performance compute cluster Euler. We describe each element of our analysis in the following.

## Energy demand projection

We project Ukraine's energy demand for the year 2060 in which Ukraine intends to be climate-neutral. We consider full electrification of the heat and transport sectors by this time as one option for complete decarbonisation. We generate an hourly time series of electricity demand in 2060 through soft-linking DESSTINEE [@Oreggioni:2022], which projects annual energy demands per sector and per energy carrier, with the Demand.ninja model [@Staffell:2023a], which upscales these demands to hourly resolution based on the prevailing weather conditions. We use the resulting profile as an exogenous input to the capacity expansion model, ignoring demand flexibility.

Annual energy demands are modelled using macroeconomic projections informing changes in energy service demands. These changes are converted to final energy consumption using assumptions on fuel switching (towards electrification) and device efficiency. We assume a 0.965% annual reduction in population between 2020--2060  [@UnitedNationsDepartmentofEconomicandSocialAffairsPopulationDivision:2022] leading to a reduction from 43.9 million people in 2020 to 29.8 million in 2060. This projection is a continuation of a thirty-year trend (Supplementary Figure\ S6). By the end of 2021 (before Russia’s invasion), Ukraine’s population had fallen 16% since its peak in 1993. This was largely due to demographics as deaths exceeded births by an average of 255,000 per year, with migration also contributing a net loss of 40,000 people per year.  Ukraine’s population fell a further 17% in 2022 due to the refugee crisis caused by Russia's invasion.  The UN expects population to rebound through to 2027, stabilising at 39 million, or 6% below the expectation from recent trends. After this point, the trend of decline is expected to continue by the UN, even in their high fertility scenario.

We assume a 4.3% annual growth in GDP per capita over the same period, taken from the high/high scenario of a report by the World Bank [@Smits:2019]. This sees GDP per capita growing from $13,000 in 2020 to $70,100 in 2060 (in Purchasing Power Parity terms). As the World Bank present a range of scenarios, we also test their low/high scenario, with 1.0% annual growth that gives $19,100 GDP per capita in 2060. Under these scenarios, Ukraine’s national GDP at PPP changes from $570\ bn in 2020 to $2,100\ bn in 2060 (high scenario) or remains flat at $570\ bn in 2060 (low scenario) as higher incomes are countered by population degrowth. We take the GDP split across economic sectors in 2020 from refs. [@WorldBankGroup:2023; @WorldBankGroup:2023a]: 70% services, 20% industry, 10% agriculture. This sectoral split is modelled to follow the trend observed across European countries, in which the service sectors grow faster than other sectors as countries get richer [@Bossmann:2015]. This yields a GDP split in 2060 of 83% services, 13% industry and 4% agriculture in the high-growth scenario, and 74% services, 18% industry and 8% agriculture in the low-growth scenario.

DESSTINEE was updated to use 2019 as the baseline calibration year, as a proxy for 2020 which we avoid due to COVID impacts on energy demand.  Final energy demands were taken from the International Energy Agency [@IEA:2023] and Statistics Service of Ukraine [@StateStatisticsServiceofUkraine:2021], shown in @tbl:final-energy-demand.  The split between end uses (for heating and cooling) is an output of the DESSTINEE model, derived from annual heating and cooling degree days.

The 2060 demand projection was constructed using a bottom-up approach in DESSTINEE. The numerous assumptions for consumer behaviour, technology mix and efficiency for each sector were derived from the default values in DESSTINEE for Ukraine’s neighbouring European countries (Poland, Slovakia, Hungary and Romania), except for assumptions on the fuel basket for each sector which were overridden to give complete electrification (in line with the study’s aims).  Influential assumptions included the average efficiency of heat pumps (Coefficient of Performance increasing from 3 in 2020 to 4.25 in 2060) and the specific energy consumption of private cars decreasing from 1.1 to 0.3\ MJ/passenger-km from 2020 to 2060 (due to efficiency gains from electric vehicles).  While these assumptions are particularly uncertain (as are any technical projections for several decades into the future), we find the results are relatively insensitive to the overall scale of demand, and thus to these particular assumptions.

In @tbl:final-energy-demand, electricity consumption is seen to increase in all sectors in both the high and low-growth scenarios.  Total energy consumption also rises in the high-growth scenario, but falls in the low-growth scenario as increasing efficiency outweighs the increase in service demands due to greater prosperity.  The residential sector is the exception to this, due to the predominance of heating demand which is modelled to see very large efficiency gains due to the uptake of heat pumps. Total demand for space heating sees a large reduction between 2020 and 2060 due to both this increasing efficiency and also the modelled increase in temperatures due to climate change: annual heating degree days are projected to fall by 17% from 3,150 to 2,600 [@Staffell:2023a]. Conversely, demand for space cooling is projected to rise sharply for two reasons.  Firstly, the existing service is already provided by air conditioners, so efficiency improvements will be marginal technology improvements, rather than fuel switching from coal and oil. Secondly, climate change increases summer temperatures and thus service demands, with annual cooling degree days growing 35% from 580 to 780 [@Staffell:2023a].

```table
---
caption: '**Historic and projected final energy demand by sector and carrier (TWh/yr).** Historic data are taken from @StateStatisticsServiceofUkraine:2021. Projections for industry do not include the production of hydrogen. {#tbl:final-energy-demand}'
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

Weather-invariant demands are modelled within DESSTINEE using the partial decomposition approach outlined in @Bossmann:2015.  Diurnal profiles (covering 24 hours) are defined for each sector and end use, for weekdays and weekends, summers and winters, and applied to the demand.  For example, heavy industry is modelled as a constant 24/7 load, residential appliances follow a general profile that mirrors household occupation.  Electric vehicles are modelled as an equal mix of time-of-use charging (daytime and evenings) and smart charging (overnight).

The Demand.ninja model was used to model weather-dependent demands for space heating and cooling in residential and commercial buildings [@Staffell:2023a].  This derives a composite temperature index from four variables which influence thermal comfort: temperature, solar irradiance, wind speeds and humidity.  From this index, daily heating and cooling degree days were derived, which are then used to apportion heating and cooling demand over time. Open-access profiles were used, derived from the model’s global average parameters, including heating and cooling thresholds of 14°C and 20°C respectively [@Demand.ninja:2023]. Demands were then upscaled from daily to hourly resolution using the diurnal profiles for heating and cooling that were derived for Ukraine’s neighbouring countries [@Staffell:2023a].

This process was first validated against historical national electricity demand in Ukraine [@ENTSO-E:2023], by applying weather data from 2018 to 2020 to the model’s baseline year (@tbl:final-energy-demand). The first ten months of 2020 were removed from comparison due to the effects of the COVID-19 pandemic on human activity and thus on energy demand [@Mehlig:2021]. This validation (@fig:load-validation) revealed a strong correlation between modelled and actual demand, with an R2\ =\ 0.95 and residual standard error of ±\ 0.5\ GW compared to mean demand of 17\ GW.

![**Estimated versus measured demand.** Demand is shown as daily averages for the years 2018--2021. In our validation, we exclude data from the COVID19 period January--October 2020 (highlighted) due to its irregular pattern.](build/results/load-validation.png){#fig:load-validation}

This same process was then applied to the 2060 scenarios, using the weather years of 2010--2014 to correspond with the years used for the wind and solar profiles (@fig:load). This profile was used in the capacity expansion model.

![**Electricity load projection used in the model.** Load is projected for the year 2060 and it is assumed that heat and transport are fully electrified. Horizontal axis shows the historical weather years used to project 2060 demand.](build/results/load.png){#fig:load}

## Capacity expansion model

To derive cost-minimal generation, storage, and transmission capacities, we apply PyPSA-Eur [@Horsch:2018] as capacity expansion model. We run the model in hourly resolution over five historical weather years (2010--2014). Given the long-term perspective in our analysis, we apply a greenfield approach in which we do not consider any existing capacities other than the transmission grid (Supplementary Figure\ S7) and hydro generation capacities. We do not consider any connections to other countries, but model Ukraine as a stand-alone system.

The model has the option to expand generation, storage, and transmission capacities of onshore wind, offshore wind, solar power, biomass, nuclear power, lithium-ion batteries, and hydrogen storage to meet demand in every hour of the year. Hydroelectricity cannot be expanded beyond pre-war’s level as the potential had been largely tapped already [@Gernaat:2017]. Fossil fuels cannot be used. Total system cost is derived by summing up annuities of investment, operation, and maintenance cost of all installed capacities (@tbl:technology-cost). While we are assessing the long-term goal of full decarbonisation by 2060, we are using cost assumptions for the year 2030 in order to provide conservative estimates for the cost of the energy system. The discount rate is set to 10% [@Andersson:2020]. The model finds the set of installed capacities with minimal total system cost.

The potential generation of solar and wind power is taken from PyPSA-Eur and has been derived in the following way.

First, the eligible area for wind and solar development is calculated for each region using *atlite* at 100\ m grid resolution [@Hofmann:2021].
For all renewable technologies, natural protection areas are excluded based on the World Database on Protected Areas (WDPA) [@UNEP-WCMC:2018].
Based on the Copernicus Global Land Cover dataset [@Buchhorn:2020], shrubland, herbaceous and sparse vegetation, and cropland are assumed to be eligible for wind and solar development.
In addition, while built-up areas are included for solar PV potentials, a distance of 1000\ m from built-up areas has to be kept for onshore wind turbines (Supplementary Figures\ S8 and S9).
Offshore wind development is allowed up to a water depth of 50\ m, which is determined based on the GEBCO bathymetry dataset [@GEBCO:2015].
Furthermore, dense shipping lanes are excluded based on the World Bank's Global Shipping Traffic Density dataset [@Cerdeiro:2021].
Wind parks further out than 30\ km from shore are assumed to be DC-connected, whereas near-shore wind parks are assumed to be AC-connected (Supplementary Figures\ S10 and S11).
For each renewable technology and region, the available area is multiplied with allowed deployment densities, approximating the socio-technical potential.
These densities are 3\ MW/km^2^ for onshore wind, 2\ MW/km^2^ for offshore wind, 1.7\ MW/km^2^ for solar.
The densities are based on purely technical estimations reported by @Scholz:2012, and our own assumption that only 30%, 20%, and 1% of eligible areas are available for onshore wind, offshore wind, and solar power.

Second, based on the eligible areas per region, capacity factor time series for wind and solar generation are calculated using *atlite*.
For this step, historical weather data for the years 2010--2014 from ECMWF's ERA5 reanalysis dataset [@Hersbach:2020] is used to convert wind speed and solar irradiance data to hourly capacity factors using models for typical wind turbines and solar panels.
The solar generation is calculated based on the incidence angle of solar irradiation, the panel tilt angle and the conversion efficiency of Crystalline Silicon (CSi) solar panels.
Power curves for a Vestas V112 3\ MW (onshore) and NREL 5\ MW (offshore) turbine are employed to map wind speeds scaled to hub height to power outputs.
The capacity factors of offshore wind generation are multiplied with a correction factor of 88.55% to approximately account for wake effects [@Bosch:2018].
Finally, the gridded dataset (0.25°x0.25°) is mapped onto the geographical shape of each region, using the available area as weighting.

We assume the future bioenergy potential to be 513\ TWh/yr (438\ TWh/yr solid biomass and 75\ TWh/yr biogas), following an assessment in @Geletukha:2020. Their projections include energy crops which may require up to ~6.6% of total land.

```table
---
caption: '**Main technology assumptions of the model using projections for 2030.** Fixed operation and maintenance cost (FOM) is given in %/year. Variable operation and maintenance cost (VOM) is given in EUR~2015~/MWh. Investment cost are given in EUR~2015~/kW. Technology lifetime is given in years. See Supplementary Note S1 for more details on the cost of nuclear power. All values are given on the AC side. Cost of solar assumes a share of 14% rooftop and 86% utility-scale solar power. {#tbl:technology-cost}'
alignment: LRRRRR
include: build/assumptions.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.25
    - 0.15
    - 0.15
    - 0.15
    - 0.15
    - 0.15
---
```

## Scenarios

We explore five scenarios to assess the impact of nuclear power on system cost of a fully-decarbonised energy supply. In the first scenario, *only-renewables*, electricity is generated exclusively from renewable sources: solar, wind, hydro, and biomass. In the second scenario, *nuclear-and-renewables*, in addition to renewables, electricity can be generated from nuclear fuel. We set the minimal level of nuclear to its pre-war level of ~55% of electricity generation: 52\ GW. In order to identify best cases for nuclear, we model nuclear capacities as being fully flexible with no ramping or minimal load constraints.

Scenarios three and four duplicate scenarios one and two but with a low-, rather than a high-economic growth assumption. To reach their pre-war level of ~55% of electricity generation, nuclear power capacities are set to 24.5\ GW or above in these scenarios.

Scenario five duplicates scenario one but with significantly lower biomass potential. While all renewable sources from scenario one are available, biomass potential is ten times smaller: 51.3\ TWh instead of 513\ TWh. Nuclear power is not available in this scenario and economic growth is high.

## Global sensitivity analysis

To understand the impact of model parameters on cost penalty of nuclear power, we perform a global sensitivity analysis. We apply the Morris Method, a method that estimates sensitivities accurately and efficiently [@Kristensen:2016; @Usher:2023]. We use the method to understand which model parameters impact the cost penalty of nuclear power the most. In the sensitivity analysis, we only consider the high-economic growth assumption scenarios. We assess ten cost parameters of the model (@tbl:gsa-parameters).

```table
---
caption: '**Parameter uncertainty ranges considered in the global sensitivity analysis.** Minimum and maximum values are relative to the default values considered in the model. {#tbl:gsa-parameters}'
alignment: LRR
include: build/gsa-parameters.csv
include-encoding: UTF-8
markdown: True
---
```

The sampling of the Morris Method works in the following way: Given ten model parameters, the method chooses a random start value within the 10-dimensional space spanned by the model parameters. The method follows a one-step-at-a-time approach, in which each iteration of the method changes the value of a single parameter only. This approach leads to a trajectory of (10 + 1) points in the parameter space. To explore this space more comprehensively, we generate 15 such trajectories, creating 165 different combinations of parameter values. We run both the scenario with and without nuclear power for all combinations allowing us to attribute cost penalty changes to parameter values. For each parameter, we observe its impact on cost penalty once for each trajectory, a total of 15 times. We use the open-source Python library SALib to perform the global sensitivity analysis [@Herman:2017].

# Results

Ukraine can supply its future decarbonised energy demand, including the shift to fully electrified heat and transport sectors, at low cost of below 90\ EUR/MWh (@fig:lcoe). This cost includes all generation, transmission, and storage capacities and does not consider trade with neighbouring countries. Requiring that nuclear capacity remains at pre-war levels incurs a cost penalty of 8% (an additional 7.2\ EUR/MWh, @fig:lcoe). To be cost competitive, nuclear capacity cost would need to be 13% lower than our assumption of about 8,000 EUR/kW (Methods). The cost and cost penalty are largely insensitive to economic growth assumptions and barely change in a low-growth scenario, in which demand is about 50% lower than the baseline (Supplementary Figure\ S1). This is despite the fact that economic growth not only impacts total energy demand, but also the shape of the load curve. Assumptions of our demand model that only alter total energy demand, such as assumptions on population count, can be expected to have an even smaller impact.

![**Cost breakdown of levelised system cost of electricity (LCOE).** Cost for scenarios with fully-electrified heat and transport sectors, average cost assumptions, and high economic growth assumption. See Supplementary Figure\ S1 for results based on both economic growth assumptions. Storage includes hydrogen, battery, and pumped-hydro storage.](build/results/lcoe-main.png){#fig:lcoe}

The largest contribution to minimal system cost stems from nuclear when it is set to pre-war levels of 55% (@fig:lcoe). Without nuclear in the energy mix, the cost of biomass is dominant (43%) while the total cost of all other renewable generation technologies plays a similar role (@fig:lcoe). Flexible balancing technologies, namely hydrogen, battery, and pumped-hydro storage, make up 17% of the system cost.

Despite its high cost contribution, biomass plays a smaller role in the electricity mixes of both minimal-cost scenarios (@fig:generation) with biomass potentials used to less than 65%. The stark discrepancy between cost and electricity contribution of biomass is caused by its role as flexibility provider in which its capacity factors are low and its specific energy cost is high. Apart from nuclear with a contribution just above 50%, solar photovoltaics generates the most electricity. Onshore wind turbines also provide significant amounts of electricity, especially without nuclear in the mix. Offshore wind turbines see almost no deployment due to their assumed high cost and sufficient onshore potentials.

![**Electricity mix.** Electricity mix of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors and with average cost assumptions and high economic growth assumption. See Supplementary Figure\ S2 for both economic growth assumptions.](build/results/generation-main.png){#fig:generation}

Meeting future energy demand with low emissions requires a substantial expansion beyond pre-war electricity generation capacities (@fig:capacities). This is especially pronounced when all generation capacities are based on renewable sources in which case more than 550\ GW need to be deployed (~250\ GW in low economic growth scenario, Supplementary Figure\ S3). Such deployment levels dwarf the ~17\ GW that existed before the start of the war and result in a required expansion of about 16\ GW/yr (7\ GW/yr) until 2060. The expansion is explained by a switch to technologies with lower capacity factors (wind and solar), growing final energy demand based on population and economic growth, and the electrification of the heat and transport sectors. Electrification and growing demand lead to a five-fold increase of pre-war electricity demand in our high-growth scenario (@fig:capacities\ D).

![**Installed capacities and electricity demand.** Installed **(A)** generation, **(B)** storage, and **(C)** energy storage capacities, and **(D)** average electricity demand of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors, average cost assumptions, and high economic growth assumption. See Supplementary Figure\ S3 for results based on both economic growth assumptions.](build/results/capacities-main.png){#fig:capacities}

In addition, the nationally-autarkic energy system we assess requires substantial amounts of energy storage capacities (@fig:capacities\ C). Storage requirements are higher when no nuclear power is in the energy mix, in which case the cost-minimal system design includes about 0.4\ TWh of battery storage capacities and almost 8\ TWh of hydrogen storage capacities in caverns. With nuclear in the energy mix, these numbers are more than halved. This can be explained by the lower amount of variable generation capacities but also by the high flexibility that we attribute to nuclear power, mitigating the need for additional flexibility from storage (Methods).

The cost penalty of keeping nuclear in the mix depends strongly on its capital cost: When capital cost is doubled, the cost penalty increases on average by more than 50\ EUR/MWh (@fig:sensitivities). The costs of biomass, solar, wind, and storage play smaller and opposite roles and decrease the cost penalty by about 10\ EUR/MWh on average for a doubling of any single model parameter (@fig:sensitivities). The impacts of these parameters strongly interact and depend on each other. For example, increasing the cost of solar can have almost no impact on cost penalty or it can have a large impact of more than -20\ EUR/MWh depending on other parameter values. This is due to the plethora of options of running a fully-renewable energy system in Ukraine with similar cost: When cost of solar increases, solar can be replaced by onshore wind, offshore wind, or biomass up to their potential limits. Similarly, if cost of biomass capacities increases, its flexibility provision can be taken over by larger hydrogen storage capacities. Only if several of these parameters are simultaneously larger than expected, there are stronger effects on the cost penalty of nuclear power.

![**Global sensitivities of cost penalty of nuclear power to model parameters based on 165 parameter combinations.** Thick blue bars show the average increase in cost penalty of keeping nuclear at least at pre-war’s levels across all parameter combinations, given a doubling of parameter values (elementary effects of the Morris method). Thin bars show the standard deviation of these effects and are therefore a proxy for the interaction between the parameters: the larger the standard deviation, the more does the impact of the parameter depend on other parameter values. The dashed grey line shows the value at which nuclear power becomes cost competitive (7.2\ EUR/MWh). All values are based on the high economic growth assumption. We ran each of the two scenarios 165 times using varying combinations of parameter values, following the Morris method.](build/results/gsa/sensitivities-annotated.png){#fig:sensitivities}

In our scenarios, both nuclear power and biomass provide significant amounts of flexibility and are able to balance fluctuations of renewable power on an annual (@fig:time-series) and diurnal (Supplementary Figures\ S4 and S5) basis. Whether nuclear power could provide such high levels of flexibility is disputed [@Morris:2018; @Jenkins:2018a] and our scenario serves therefore as a best case for nuclear: Should nuclear capacities be unable to provide these required levels of flexibility, other flexibility  providers had to be included for extra cost. The role of biomass is less unique: If less biomass is available than we assume, its flexibility provision can be replaced by storage capacities (right panel in @fig:time-series) at low extra cost (@fig:sensitivities).

![**Generation time series across all assessed weather years 2010--2014.** The areas show monthly generation of the different components of the energy system. See Supplementary Figures\ S4 and S5 for zooms into selected summer and winter weeks. In the low-biomass scenario, available biomass is constrained to 10% of the otherwise assumed value of 513\ TWh/yr. Storage includes hydrogen, battery, and pumped-hydro storage. Storage demand is not shown.](build/results/time-series-main.png){#fig:time-series}

# Discussion

We show that a decarbonised energy supply in Ukraine is technically feasible and economically viable. With electrified heat and transport sectors, the combination of renewable generation and energy storage technologies can supply electricity at a system cost of below 90\ EUR/MWh, including generation, storage, and transmission. The projected transition to an energy system without fossil fuels can fully meet the growing demand and full electrification of industrial and residential energy consumption, including mobility, according to various forecasts of economic growth for Ukraine. Building such a system that is able to cope with population, economic, and demand growth, however, requires a substantial and sustained effort deploying new infrastructure.

We further show that none of the technologies is essential on its own, leaving room for decision making based on risks, uncertainties, and societal or political preferences. Nuclear power is not essential and only appears in the cost-minimal energy mix when we prescribe its use. Based on today’s cost estimate (7,940\ EUR/kW), it increases rather than decreases cost of energy [@Steigerwald:2023]. Solar power takes on a prominent role in the cost-minimal energy mixes due to its low and falling cost, but, as we show, it can be replaced by more wind or biomass should its cost developments be less favourable. Similarly, our cost-minimal energy mixes contain large amounts of electricity from biomass power stations, but their flexibility provision can be replaced by larger energy storage capacity. There is not just one option, but a plethora of options to build a decarbonised energy system for Ukraine, similar as has been shown for other parts of the world [@Pickering:2022].

The potential of renewable energy in Ukraine far exceeds the estimated demand by 2060: even assuming high economic growth, the potentials of solar and wind power are roughly three times larger than demand. This indicates that Ukraine has the potential to sell surplus energy generated. By using opportunities tendered by the newly established Ukraine-EU electricity interconnection [@ENTSO-E:2022a], Ukraine can transit from being a net importer of energy resources to being a net exporter in the post-war time. Ukraine's geographic location places it at the crossroads of major energy transit routes, including pipelines and transportation networks. This strategic position allows Ukraine to serve as a transit country for energy resources moving between Europe and Asia. More extensive trade, in turn, will likely reduce the potential cost of electricity for domestic consumers.

Our analysis was deliberately conservative with regards to a fully renewable electricity supply and could overestimate flexibility needs. We therefore expect its levelised cost to be an upper bound, and the levelised cost penalty to be a lower bound for nuclear power with respect to flexibility needs. This is for three reasons.

First, we model Ukraine’s electricity supply in isolation. As has been shown in other studies [@Trondle:2020a; @Neumann:2021; @Brown:2021a], small-scale, autarkic systems show higher flexibility needs. By no means is our model choice based on a political suggestion, but it serves as an extreme case in which both flexibility needs and cost are higher than in a case in which the electricity system is connected to its neighbours. Using this model choice, we can safely assume that the results of both scenarios are pessimistic in this regard.

Second, our analysis ignores flexibility provision from outside the electricity system therefore increasing flexibility needs within it. Even in a fully-electrified scenario as ours, other sectors could provide flexibility, e.g. demand-side responses for heat storage or battery-electric mobility. Here as well, our analysis yields an upper bound in this regard.

Third, our analysis assumes full flexibility of nuclear power, thereby decreasing flexibility needs from other sources in the nuclear scenario. In our model, nuclear power has no constraints on ramping behaviour or on minimal loads and is therefore able to follow and balance fluctuations from renewable sources perfectly. As this is not possible with current technology, our approach yields a lower bound for the cost penalty of nuclear power.

Building such a decarbonised energy system requires a significant infrastructural effort: depending on the assumption on economic growth and whether nuclear is part of the mix, it requires more than 120--550\ GW of newly deployed renewable generation infrastructure, 23--80\ GW of electricity storage infrastructure, and 0--50\ GW of nuclear capacities. To reach a carbon-neutral electricity system by 2060, this requires an average annual deployment of 3.5--16\ GW of renewable generation infrastructure which is far above most recent national goals of about 1\ GW per year and unlikely to occur without long-lasting financial and regulatory political support.

Apart from financial and regulatory support, coordination of this infrastructural effort will be key. Not only large parts of the generation infrastructure from the time before Russia’s invasion are destroyed, but also large parts of the transmission infrastructure. The parallel build-up of both generators and grid will require thoughtful coordination. This is not only a challenge but also an opportunity to develop a modern electricity grid that is designed to handle large amounts of distributed and variable energy sources. Possibly, the build-up could be orchestrated from the bottom to the top, focussing on restoring supply on the local level and leading to a supply that is stronger decentralised [@Bauknecht:2020] and relying more on bioenergy [@Mangoyana:2011]. Whether more or less decentralised, such a modern energy supply could provide a blueprint for decarbonised energy provision.

Our study has implications for Ukraine, the European Union, and beyond. Rebuilding Ukraine’s energy supply in a decarbonised way will help reduce greenhouse gas emissions, improve environmental quality, create jobs [@Trypolska:2021], and increase energy security. The open data and models we share enables possibilities for broad research tackling these and further aspects of the energy system of Ukraine. For the European Union, coordinated efforts with Ukraine may contribute to the achievement of climate goals. Beyond Europe, building on the literature on green recovery of economic crises [@Bersalli:2023], our study showcases the potential that windows of opportunity may provide globally.

# Lead Contact

Requests for further information, resources, and materials should be directed to the lead contact, Tim Tröndle (tim.troendle@usys.ethz.ch).

# Data and Code Availability

The datasets generated during this study are available on Zenodo with DOI\ [10.5281/zenodo.10814894](https://doi.org/10.5281/zenodo.10814894).

We also refer to the documentation of the PyPSA energy system modelling framework ([https://pypsa.readthedocs.io](https://pypsa.readthedocs.io)) and of the European energy system model PyPSA-Eur ([pypsa-eur.readthedocs.io](https://pypsa-eur.readthedocs.io)) for technical instructions on how to
install, modify and run the model as well as more detailed explanations of the model structure.

The code of the specific version used, the configuration to reproduce the results and all analysis steps are publicly available as a reproducible Snakemake [@Koster:2012] workflow on Zenodo with DOI [10.5281/zenodo.11276962](https://doi.org/10.5281/zenodo.11276962).

# CRediT author statement
**Conceptualisation**, T.T., and A.P.; **Methodology**, T.T., F.N., and I.S.; **Investigation**, T.T.; **Data curation**, T.T., O.M, O.T., F.N., and I.S.; **Formal analysis**, T.T., F.N., I.S.; **Software**, T.T., F.N., I.S.; **Visualisation**, T.T., and F.N.; **Writing -- original draft**, T.T., O.M, F.N., and I.S.; **Writing -- review & editing**, O.T., V.P., and A.P.; **Funding acquisition**: A.P.

# Declaration of interests
The authors declare no competing interests.

# Acknowledgements
This work has been supported by a Scholars at Risk Grant by the Swiss National Science Foundation (SNSF).

# References
