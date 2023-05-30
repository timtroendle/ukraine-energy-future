# Introduction

We are analysing different options for Ukraine's future electricity supply.

# Methods

I am using PyPSA-Eur [@Horsch:2018] to model Ukraine's electrictiy supply. Its main assumptions about individual technologies are shown in @tbl:assumptions. The depreciation rate is set to 7%.

For now, I am using historic electricity demand from 2013.

For now, I run two scenarios. In the first scenario, `only-renewables`, electricity is generated exclusively from renewables: solar, wind, hydro, and biomass. The lion's share of generation stems from solar and wind.

In the second scenario, `nuclear-and-renewables`, in addition to renewables, electricity can be generated from ... coal. It should be nuclear really, but for now the model has an issue with that and so I included coal instead of nuclear. Coal generation has similar technical characteristics but different costs and I will fix this as soon as possible.

```table
---
caption: 'Main technology assumptions of the model. Fixed operation and maintenance cost (FOM) is given in "% / year". Variable operation and maintenance cost (VOM) is given in "EUR / MWh". Investment cost are given in "EUR / kW". Technology lifetime is given in years. {#tbl:assumptions}'
alignment: LRRRR
include: build/assumptions.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.1
    - 0.1
    - 0.1
    - 0.1
---
```

# Results

@tbl:lcoe shows the levelised cost of electricity (including cost of the transmission system), @tbl:capacities-power shows the installed generation capacities in the scenarios, and @tbl:capacities-energy shows the installed storage capacities in the scenarios.

```table
---
caption: 'Levelised cost of electricity. {#tbl:lcoe}'
alignment: LR
include: build/results/lcoe.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.3
    - 0.2
---
```

```table
---
caption: 'Installed generation capacities (GW). {#tbl:capacities-power}'
alignment: LRR
include: build/results/capacities-power.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.3
    - 0.3
---
```

```table
---
caption: 'Installed storage capacities (GWh). {#tbl:capacities-energy}'
alignment: LRR
include: build/results/capacities-energy.csv
include-encoding: UTF-8
markdown: True
width:
    - 0.2
    - 0.3
    - 0.3
---
```

# References
