####################################################################
# Author: Franny Dean
# Purpose: Model incorporating assumption that P_c decays linearly
# Date: 4/9/2021
#
####################################################################

from typing import Dict
from numbers import Number
from functools import partial
import numpy as np
from scipy.optimize import bisect


class MIRModel:
    """Mortality incidence ratio model to approximate the survival rate. This model
       assumes that i, other_mortality are constant over time, and that probability of 
       death due to cancer decays exponentially. 
    """

    def __init__(self,
                 mir: np.ndarray,
                 other_mortality: np.ndarray,
                 disease_period: np.ndarray,
                 slope: np.ndarray):
        """Constructor of MIRModel.

        Args:
            mir (np.ndarray): Mortality and incidence ratio.
            other_mortality (np.ndarray): Probablity of death from other causes.
            disease_period (np.ndarray):
                Excess mortality period. Either to be a positive integer or
                ``np.inf``. Expected to be cause specific.
            slope (np.ndarray): the slope associated to the exponential linear decay with time 
                of P_c. Expected to be cause specific.
        """
        # parse inputs
        self.mir = np.asarray(mir)
        self.other_mortality = np.asarray(other_mortality)
        if np.isinf(disease_period):
            self.disease_period = disease_period
        else:
            self.disease_period = int(disease_period)
        self.num_points = self.mir.size
        self.slope = np.asarray(slope)

        # assertions for the inputs
        self._check_inputs()

        # place holder for base excess mortality
        self.base_excess_mortality = np.full(self.num_points, np.nan)

    def _check_inputs(self):
        """Check the inputs for the class.
        """
        assert self.mir.size == self.num_points
        assert ((self.mir >= 0) & (self.mir <= 1)).all()

        assert self.other_mortality.size == self.num_points
        assert (self.other_mortality >= 0).all()
        assert (self.other_mortality <= 1).all()

        assert self.disease_period.size == self.num_points
        assert (self.disease_period >= 1).all()

        assert self.slope.size == self.num_points

    def compute_base_excess_mortality(self):
        """Compute the base probability of death from the disease in year of diagnosis (P_c(0)).
        """

        for i in range(self.num_points):
            if self.mir[i] >= 1 - self.other_mortality[i]:
                self.base_excess_mortality[i] = 1 - self.other_mortality[i]
            else:
                self.base_excess_mortality[i] = bisect(
                    partial(self.mir_equation,
                            other_mortality=self.other_mortality[i],
                            mir=self.mir[i],
                            disease_period=self.disease_period,
                            slope=self.slope[i]),
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

        if np.isnan(self.base_excess_mortality).all():
            print("Compute base excess mortality first.")
            self.compute_base_excess_mortality()

        abs_survival_rate = cumulative_survival(self.base_excess_mortality, self.slope, num_years, self.other_mortality)
        ref_survival_rate = (1 - self.other_mortality)**num_years
        rel_survival_rate = abs_survival_rate/ref_survival_rate

        # Would it be better to calculate the cause specific survival? 1 - P_c(t)?

        return {'abs': abs_survival_rate, 'rel': rel_survival_rate}

    @staticmethod
    def mir_equation(base_excess_mortality: float,
                     other_mortality: float,
                     mir: float,
                     disease_period: int, 
                     slope: float) -> float:
        """Equation that use MI ratio approximating the excess mortality (P_c).

        Args:
            base_excess_mortality (float): Probablity of death from the disease in year of diagnosis.
            other_mortality (float): Probablity of death from other causes.
            mir (float): Mortality incidence ratio.
            disease_period (int): Excess mortality period.
            slope (float): the cause-specific slope associated to the exponential linear decay with time of P_c.

        Returns:
            float: Value of the difference between the expected and the truth.
        """
        val_true = 0
        for year in list(range(1,disease_period)):
            cumulative_survival = cumulative_survival(base_excess_mortality, slope, year, other_mortality)
            val_true += P_c(base_excess_mortality, slope, year)*cumulative_survival
        
        
        val_expected = mir

        return val_true - val_expected
    
    @staticmethod
    def P_c(base_excess_mortality: float,
            slope: float, 
            year: int) -> float:
        """Linear decay of the probability of death.

        Args:
            base_excess_mortality (float): Probablity of death from the disease in year of 
                diagnosis.
            slope (float): the cause-specific slope associated to the exponential linear decay 
                with time of P_c.
            year (int): the whole years since diagnosis.

        Returns:
            float: Returns the year specific excess mortality.
        """

        value = base_excess_mortality*exp(slope*year)

        return value
    
    @staticmethod
    def cumulative_survival(base_excess_mortality: float,
            slope: float, 
            year: int,
            other_mortality: float) -> float:
        """Calculates cumulative survival based on probability of death.

        Args:
            base_excess_mortality (float): Probablity of death from the disease in year of 
                diagnosis.
            slope (float): the cause-specific slope associated to the linear exponential decay with time 
                of P_c.
            year (int): the whole years since diagnosis.
            other_mortality (float): Probablity of death from other causes.

        Returns:
            float: Returns the year specific cumulative survival.
        """
        value = 1
        for i in list(range(1,year)):
            value *= 1 - P_c(base_excess_mortality, slope, i) - background_mortality
        
        

        return value

