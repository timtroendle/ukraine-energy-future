OPTIMALITY_CONDITION = "INFO:pypsa.linopf:Optimization successful."
SUBOPTIMALITY_CONDITION = "WARNING:pypsa.linopf:Optimization solution is sub-optimal."


def test_scenario_is_optimal(scenario_log: list[str]):
    assert any(OPTIMALITY_CONDITION in line for line in scenario_log)


def test_scenario_is_not_suboptimal(scenario_log: list[str]):
    assert not any(SUBOPTIMALITY_CONDITION in line for line in scenario_log)
