import pandas as pd


TOO_LARGE_CAPACITY_FACTOR = 0.6
TOO_LOW_CAPACITY_FACTOR = 0.09


def test_capacity_factors_not_too_large(re_potential_time_series: pd.Series):
    avg = re_potential_time_series.mean().item()
    assert avg < TOO_LARGE_CAPACITY_FACTOR


def test_capacity_factors_not_too_low(re_potential_time_series: pd.Series):
    avg = re_potential_time_series.mean().item()
    assert avg > TOO_LOW_CAPACITY_FACTOR
