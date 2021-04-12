"""
    model
    ~~~~~
"""
from typing import Dict
from numbers import Number
from functools import partial
import numpy as np
from scipy.optimize import bisect


class MIRModel:
    """Mortality incidence ratio model to approximate the survival rate.
    """

    def __init__(self,
                 mir: np.ndarray,
                 other_mortality: np.ndarray,
                 disease_period: float = 5):
        """Constructor of MIRModel.

        Args:
            mir (np.ndarray): Mortality and incidence ratio.
            other_mortality (np.ndarray): Probablity of death from other causes.
            disease_period (float, optional):
                Excess mortality period. Either to be a positive integer or
                ``np.inf``. Defaults to 5.
        """
        # parse inputs
        self.mir = np.asarray(mir)
        self.other_mortality = np.asarray(other_mortality)
        if np.isinf(disease_period):
            self.disease_period = disease_period
        else:
            self.disease_period = int(disease_period)
        self.num_points = self.mir.size

        # assertions for the inputs
        self._check_inputs()

        # place holder for excess mortality
        self.excess_mortality = np.full(self.num_points, np.nan)

    def _check_inputs(self):
        """Check the inputs for the class.
        """
        assert self.mir.size == self.num_points
        assert ((self.mir >= 0) & (self.mir <= 1)).all()

        assert self.other_mortality.size == self.num_points
        assert (self.other_mortality >= 0).all()
        assert (self.other_mortality <= 1).all()

        assert isinstance(self.disease_period, Number)
        assert self.disease_period >= 1

    def compute_excess_mortality(self):
        """Compute the probability of death from the disease.
        """
        if np.isinf(self.disease_period):
            excess_mortality = self.other_mortality * \
                self.mir/(1 - self.mir)
            self.excess_mortality = np.maximum(
                0, np.minimum(1 - self.other_mortality, excess_mortality))

        else:
            for i in range(self.num_points):
                if self.mir[i] >= 1 - self.other_mortality[i]:
                    self.excess_mortality[i] = 1 - self.other_mortality[i]
                else:
                    self.excess_mortality[i] = bisect(
                        partial(self.mir_equation,
                                other_mortality=self.other_mortality[i],
                                mir=self.mir[i],
                                disease_period=self.disease_period),
                        0.0, 1.0
                    )

    def get_survival_rate(self,
                          num_years: int) -> Dict[str, np.ndarray]:
        """Absolute survival rate.

        Args:
            num_years: Number of years for the survival rate.

        Returns:
            Dict[str, np.ndarray]:
                Calculate the absolute and relative survival rate from MI ratio,
                and return them in a dictionary, with key ``'abs'`` for the
                absolute and key ``'rel'`` for the relative.
        """
        if num_years < 1:
            raise ValueError(f"Number of years of survival must be greater or"
                             f"equal than 1, was {num_years}.")

        if np.isnan(self.excess_mortality).all():
            print("Compute excess mortality first.")
            self.compute_excess_mortality()

        prob_death_per_year = self.excess_mortality + self.other_mortality
        abs_survival_rate = (1 - prob_death_per_year)**num_years
        ref_survival_rate = (1 - self.other_mortality)**num_years
        rel_survival_rate = abs_survival_rate/ref_survival_rate

        return {'abs': abs_survival_rate, 'rel': rel_survival_rate}

    @staticmethod
    def mir_equation(excess_mortality: float,
                     other_mortality: float,
                     mir: float,
                     disease_period: int) -> float:
        """Equation that use MI ratio approximating the excess mortality.

        Args:
            excess_mortality (float): Probablity of death from the disease.
            other_mortality (float): Probablity of death from other causes.
            mir (float): Mortality incidence ratio.
            disease_period (int): Excess mortality period.

        Returns:
            float: Value of the difference between the expected and the truth.
        """
        prob_death = excess_mortality + other_mortality
        prob_survival = 1 - prob_death

        val_true = excess_mortality*(1 - prob_survival**disease_period)
        val_expected = mir*prob_death

        return val_true - val_expected
