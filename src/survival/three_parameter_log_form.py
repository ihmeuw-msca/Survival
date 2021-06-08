####################################################################
# Author: Franny Dean
# Purpose: Model incorporating assumption that log(P_c) follows b*log(1+exp(a-x))+c
# Date: June 3, 2021
#
####################################################################

from typing import Dict
from numbers import Number
from functools import partial
import numpy as np
from scipy.optimize import bisect
import math
import scipy.special as scipy


def cumulative_survival(a: float,
            b: float,
            c: float,
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
        for i in range(1,year):
            
            value *= (1 - P_c(a, b, c, i) - other_mortality)
        
        return value


def P_c(a: float,
        b: float,
        c: float, 
        year: int) -> float:
        """Linear decay of the probability of death.

        Args:
            a, b, c (floats): the parameters associated with the model b*log(1-exp(x-a))+c
            year (int): the whole years since diagnosis.

        Returns:
            float: Returns the year specific excess mortality.
        """
        value = scipy.expit(b*math.log(1+math.exp(-(year-a)))+c)
        return value


class LogFormMIRModel:
    """Mortality incidence ratio model to approximate the survival rate. This model
       assumes that i, other_mortality are constant over time, and that the log of probability of 
       death due to cancer changes according to b*log(1+exp(a-x))+c for some a,b,c. 
    """

    def __init__(self,
                 mir: np.ndarray,
                 other_mortality: np.ndarray,
                 a: np.ndarray,
                 b: np.ndarray):
        """Constructor of MIRModel.

        Args:
            mir (np.ndarray): Mortality and incidence ratio.
            other_mortality (np.ndarray): Probablity of death from other causes.
            a, b (np.ndarrays): parameters associated with the functional form for log probability of death
        """
        # parse inputs
        self.mir = np.asarray(mir)
        self.other_mortality = np.asarray(other_mortality)
        self.a = np.asarray(a)
        self.b = np.asarray(b)
        self.num_points = self.mir.size
        
        # assertions for the inputs
        self._check_inputs()

        # place holder for c parameter 
        self.c = np.full(self.num_points, np.nan)
        

       
    def _check_inputs(self):
        """Check the inputs for the class.
        """
        assert self.mir.size == self.num_points
        assert ((self.mir >= 0) & (self.mir <= 1)).all()

        assert self.other_mortality.size == self.num_points
        assert (self.other_mortality >= 0).all()
        assert (self.other_mortality <= 1).all()

        assert self.a.size == self.num_points
        assert self.b.size == self.num_points
        
    def compute_third_parameter(self):
        """Compute the base probability of death from the disease in year of diagnosis (P_c(0)).
        """

        for i in range(self.num_points):
            
            if self.mir[i] >= 1 - self.other_mortality[i]:
                # this might be wrong?
                self.c[i] = scipy.logit(1 - self.other_mortality[i])-self.b[i]*math.log(1+math.exp(self.a[i]))
            else:
                self.c[i] = bisect(    
                    partial(self.mir_equation,
                            other_mortality=self.other_mortality[i],
                            mir=self.mir[i],
                            a=self.a[i],
                            b=self.b[i]),
                    -10.0, 0.0
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

        if np.isnan(self.c).all():
            print("Compute third parameter first.")
            self.compute_third_parameter()

        abs_survival_rate = [1]*self.num_points
        rel_survival_rate = [1]*self.num_points
        for i in range(self.num_points):
            abs_survival_rate[i] = cumulative_survival(self.a[i], 
                                                       self.b[i], 
                                                       self.c[i],
                                                       num_years, 
                                                       self.other_mortality[i])
            rel_survival_rate[i] = abs_survival_rate[i]/((1-self.other_mortality[i])**num_years)
            # "Net" survival.
            #for t in range(1, num_years):
            #    rel_survival_rate[i] *= (1 - P_c(self.a[i], self.b[i], self.c[i], t))
        
        
        return {'abs': abs_survival_rate, 'rel': rel_survival_rate}
    
    
    
    @staticmethod
    def mir_equation(c: float,
                     other_mortality: float,
                     mir: float,
                     a: float,
                     b: float) -> float:
        """Equation that use MI ratio approximating the excess mortality (P_c).

        Args:
            other_mortality (float): Probablity of death from other causes.
            mir (float): Mortality incidence ratio.
            a, b, c (floats): parameters associated to the model
        Returns:
            float: Value of the difference between the expected and the truth.
        """
        val_true = 0
        
        
        for year in range(1, 16):
            cum_survival = cumulative_survival(a, b, c, year, other_mortality)
            val_true += P_c(a, b, c, year)*cum_survival
     
        val_expected = mir
        
        return val_true - val_expected
    
  
 
  