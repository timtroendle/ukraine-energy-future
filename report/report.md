# Introduction

We are analysing different options for Ukraine's future electricity supply.

# Methods

I am using PyPSA-Eur to model Ukraine [@Horsch:2018].

For now, I am using historic electricity demand from 2013.

For now, I run two scenarios. In the first scenario, `only-renewables`, electricity is generated exclusively from renewables: solar, wind, hydro, and biomass. The lion's share of generation stems from solar and wind.

In the second scenario, `nuclear-and-renewables`, in addition to renewables, electricity can be generated from ... coal. It should be nuclear really, but for now the model has an issue with that and so I included coal instead of nucelar. Coal generation has similar technical characteristics (albeit different costs) but I will fix this as soon as possible.

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
    - 0.3
    - 0.2
    - 0.2
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
    - 0.3
    - 0.2
    - 0.2
---
```

# References
