import pytest


OPTIMALITY_CONDITION = "INFO:pypsa.linopf:Optimization successful."
SUBOPTIMALITY_CONDITION = "WARNING:pypsa.linopf:Optimization solution is sub-optimal."


def test_main_scenario_is_optimal(main_scenario_log: list[str]):
    assert any(OPTIMALITY_CONDITION in line for line in main_scenario_log)


@pytest.mark.xfail # FIXME ensure gsa runs are not suboptimal
def test_gsa_scenario_is_optimal(gsa_scenario_log: list[str]):
    assert any(OPTIMALITY_CONDITION in line for line in gsa_scenario_log)


def test_main_scenario_is_not_suboptimal(main_scenario_log: list[str]):
    assert not any(SUBOPTIMALITY_CONDITION in line for line in main_scenario_log)


@pytest.mark.xfail # FIXME ensure gsa runs are not suboptimal
def test_gsa_scenario_is_not_suboptimal(gsa_scenario_log: list[str]):
    assert not any(SUBOPTIMALITY_CONDITION in line for line in gsa_scenario_log)
