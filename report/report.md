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

System cost of renewable only scenario with fully-electrified heat and transport sectors and with average cost and medium-high economic growth assumption is low (<70\ €/MWh, @tbl:lcoe). Keeping nuclear at pre-war's levels shows 18% (12.2\ €/MWh) higher cost (@tbl:lcoe). System cost includes generation, transmission, and storage capacities.

```table
---
caption: 'Levelised cost of electricity. {#tbl:lcoe}'
alignment: LR
include: build/results/lcoe.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.4
    - 0.2
---
```

Meeting future electricity demand requires a substantial expansion of pre-war generation capacities in both scenarios (@fig:capacities-power), but especially in the renewables only scenario. This is explained by a switch to technologies with lower capacity factors (wind and solar), growing final energy demand based on population and economic growth, and an electrification of the heat and transport sectors (==can we compare final energy demand pre-war and scenarios?==).

![**Installed capacities.** Installed **(A)** generation and **(B)** storage capacities of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors and with average cost assumptions and medium-high economic growth assumption.](build/results/capacities-power.png){#fig:capacities-power}

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

![**Electricity mix.** Electricity mix of the pre-war electricity system and for scenarios with fully-electrified heat and transport sectors and with average cost assumptions and medium-high economic growth assumption.](build/results/generation.png){#fig:generation}

The cost of nuclear has the greatest impact on LCOE difference between the two scenarios (@fig:sensitivities). When Nuclear capital cost is doubled, the nuclear-and-renewables scenario becomes on average more than 90\ €/MWh more expensive compared with the only-renewables scenario. Doubling the capital cost of solar, wind, and biomass increases the relative cost of the only-renewables scenario, but to a lower degree (20--30\ €/MWh on average). Capital costs of storage technologies and marginal costs have low or almost no impact.

![**Sensitivities of LCOE difference between scenarios to model parameters.** Thick bars show the average increase in LCOE of the nuclear-and-renewables compared with the only-renewables scenario, given a doubling of parameter values. Thin bars show the standard deviation of these effects and are therefore a proxy for the interactivity of the parameters: the larger the standard deviation, the more does the impact of the parameter depend on other parameter values.](build/results/gsa/sensitivities.png){#fig:sensitivities}

# Discussion

Both scenarios require large amounts of new generation and storage capacities. Given that biomass fuel is currently not constrained in the model, we can assume that storage capacities will increase even further once biomass fuel is constrained. However, we can also assume that storage capacities can decrease strongly with connecting Ukraine to neighbouring countries.

While required generation and storage capacities are large, there is no reason to believe these capacities could not be built. Ukraine can be powered by only renewables in the future at cost similar to a energy supply including nuclear.

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

We project Ukraine's energy demand for the year 2060 in which Ukraine intends to be climate-neutral. We consider full electrification of the heat and transport sectors by this time as one option for complete decarbonisation. We generate an hourly time series for a single weather year (2013 ==Iain, why this year?==) which we use as a exogenous input for the capacity expansion model, ignoring demand flexibility.

We project annual energy demand in 2060 using DESSTINEE [@Oreggioni:2022] which projects energy demand per sector and per energy carrier. For the growth of the individual sectors, we apply the following assumptions. We assume a ==ZZ%== annual population growth between 2020--2060 (==SOURCE==) leading to a population of ==XX== million people  in 2060. We assume a ==ZZ%== annual GDP per capita growth in the same duration (==SOURCE==), leading to GDP per capita growing from ==XX== to ==YY==. We assume the GDP split across the different sectors to change in a way in which services grows fastest as countries get richer, leading to a GDP split in 2060 of ==XX%== services, ==YY%== agriculture, and ==ZZ%== industry (==SOURCE==).

Using these sector growth projections and historical final energy demand (@tbl:final-energy-demand-2020), DESSTINEE projects future final energy demand per sector. This includes a full electrification and energy efficiency improvements. The projected electricity demand in 2060 sums to ==YY== TWh/yr.

-----------------------------------------------------------
Sector                            Electricity      Other fuels
---------------      ------------------------ ----------------
Industry                      46                   140

Commercial                   22                   35

Residential                   37                   121

Transport                     6                    87

Other                         4                    16
(agriculture,
fishing)
------------------------------------------------------------

: Final energy demand by sector and carrier (TWh/yr in 2020). ==Olena T.: SOURCE?== {#tbl:final-energy-demand-2020}

To derive a hourly load time series of the fully-electrified energy demand, we disaggregate the annual energy demand in time. For commercial and residential heat, we use demand.ninja ==(CITE ONCE PUBLISHED)==. For transport we do X ==Iain: what?==. For electricity, we do Y ==Iain: what?==. For industry, we do Z ==Iain: what?==. Summing all these elements leads to a national load time series (@fig:load) which we feed into the capacity expansion model.

![**Electricity load projection used in the model (preliminary).** Load is project for the year 2060 and it is assumed that heat and transport are fully electrified.](build/results/load.png){#fig:load}

## Capacity expansion model

To derive cost-minimal generation, storage, and transmission capacities, we apply PyPSA-Eur [@Horsch:2018] as capacity expansion model. Given the long-term perspective in our analysis, we apply a greenfield approach in which we do not consider any existing capacities other than the transmission grid and hydro generation capacities. We do not consider any connections to other countries, but model Ukraine as a stand-alone system.

The model has the option to expand generation and storage capacities of onshore wind, solar power, biomass, nuclear power, lithium-ion batteries, and hydrogen storage to meet demand in every hour of the year. Total system cost is derived by summing up annuities of investment, operation, and maintenance cost of all installed capacities (@tbl:technology-cost). The depreciation rate is set to 10% [@Andersson:2020]. The model finds the set of installed capacities with minimal total system cost.

The potential generation of solar and wind power is taken from PyPSA-Eur and has been derived in the following way. ==Fabian, can you explain?==

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
