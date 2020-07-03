"""
    test_model
    ~~~~~~~~~~
"""
import numpy as np
import pytest
from survival import MIRModel


num_points = 5
deaths = np.random.rand(num_points)*100.0 + 1.0
cases = deaths + np.random.rand(num_points)*1000.0 + 100.0
mir = deaths/cases
other_mortality = np.random.rand(num_points)*0.05


@pytest.fixture
def model():
    return MIRModel(mir, other_mortality, disease_period=5)


def test_compute_excess_mortality(model):
    model.compute_excess_mortality()
    assert np.allclose([
        model.mir_equation(model.excess_mortality[i],
                           model.other_mortality[i],
                           model.mir[i],
                           model.disease_period)
        for i in range(model.num_points)
    ], 0.0)


def test_get_survival_rate(model):
    model.compute_excess_mortality()
    result = model.get_survival_rate(num_years=5)
    assert result['abs'].size == model.num_points
    assert result['rel'].size == model.num_points
