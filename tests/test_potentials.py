def test_biomass_generation_is_positive(biomass_generation: float):
    assert biomass_generation >= 0


def test_biomass_generation_below_potential(biomass_potential: float, biomass_generation: float):
    assert biomass_generation <= biomass_potential
