##################################################
# Author: Franny Dean
# Purpose: Age group extension to survival model
# Date: 9/1/2020
#
##################################################

from typing import Dict
from numbers import Number
from functools import partial
import pandas as pd
import numpy as np 
from scipy.optimize import bisect

class AgeSurvivalModel:
    """ Mortality incidence ratio used to approximate age specific survival rate
    """

    def _init_(self, inputs: pd.DataFrame):
        """Constructor of the model.

            Arguments:
                age_id (np.ndarray): the list of the corresponding age groups
                mir (np.ndarray): mortality to incidence ratios for an entire age series
                other_mortality (np.ndarray): probability of death from other causes for an entire age series
        """
        required = ['age_group_id','location_id','other_mortality','sex_id','acause','year_id','mi_ratio']
        assert [column for column in required if(ele in self.inputs.columns.tolist())]
        #HELP: I think something the syntax in the above is wrong?

        self.inputs = pd.DataFrame(inputs[required])
        self.num_points = len(inputs.axes[0])
        
        #check inputs variable values
        self._check_inputs()

        #place holders
        self.inputs['P_c'] = np.full(self.num_points, np.nan)
        self.inputs['P_s'] = np.full(self.num_points, np.nan)
        self.inputs['abs_survival'] = np.full(self.num_points, np.nan)
        self.inputs['rel_survival'] = np.full(self.num_points, np.nan)

        #create UID column on location-year-sex-cause
        self.inputs['UID'] = str(self.inputs['location_id']+
                                str(self.inputs['year_id'])+
                                str(self.inputs['sex_id'])+
                                self.inputs['acause'])

        #create df dictionaries
        self.input_dfs = dict()
        
        for uid in self.inputs['UID'].unique():
            uid_data = self.inputs[inputs['UID'] == uid].copy()
            #sort by age group (ascending)
            uid_data.sort_values(by=['age_group_id'], inplace=True)
            self.input_dfs[uid] = uid_data


    def _check_inputs(self):
        """Check the values of the inputs for the class for suitability.
        """
        #no NAs
        self.inputs.dropna()

        #want other_mortality less than 1, greater than zero
        assert (self.inputs['other_mortality']>=0).all()
        assert (self.inputs['other_mortality']<1).all()

    def compute_P_s(self):
        """Main function which computes the probablity of survival P_s for each
            age group in sequence.
        """
        for uid, data in self.input_dfs.items():
            #compute base case age groups
            self.calculate_survival_first_age_group(data, self.data.age_group_id.iloc[0])
            self.calculate_survival_second_age_group(data, self.data.age_group_id.iloc[1],self.data.age_group_id.iloc[0])

            #compute remaining age groups
            for i in range(2,age_id.size()-1):
                calculate_survival_other_age_group(data,
                                                    self.data.age_group_id.iloc[i],
                                                    self.data.age_group_id.iloc[i-1],
                                                    self.data.age_group_id.iloc[i-2])

    def compute_n_year_survival(self, num_years: int):
        """Main function for getting the n-year relative and absolute survival.
        """
        #validate input years of survival
        if num_years < 1:
            raise ValueError(f"Number of years of survival must be greater or"
                             f"equal than 1, was {num_years}.")

        #catch to make sure P_s has been calculated
        if np.isnan(self.inputs['P_s']).all():
            print("Computing survival first.")
            self.compute_P_s()

        #calculate survival directly from P_s
        for i in range(1,self.num_points):
            self.inputs.abs_survival[i] = self.inputs.P_s[i]**num_years
            self.inputs.rel_survival[i] = self.inputs.abs_survival[i]/(1-self.inputs['other_mortality'][i]**num_years)

    ### SOLVE EQUATIONS ###
    def calculate_survival_first_age_group(self, uid_data: pd.DataFrame, age_id: int):
        """The youngest age group (below which there is assumed to be no incidence).
        """
        #find first age group position
        i = np.where(uid_data.age_group_id==age_id)
        
        #check to make sure that MIR is within bounds
        if uid_data.mi_ratio[i] < (1-uid_data.other_mortality[i]):

            #solve for P_c (probability of death due to cause)
            uid_data.P_c[i] = bisect(
                    partial(self.first_age_group_equation, 
                    other_mortality=uid_data.other_mortality[i],
                    mir=uid_data.mir[i]), 0.0, 1.0)
        else:
            #MIR is too large to properly solve (survival will go negative) assign survival of 0
            uid_data.P_c[i] = 1 - uid_data.other_mortality


        #solve for P_s (probability of abs survival in 1 year)
        uid_data.P_s[i] = 1 - (uid_data.P_c[i] + uid_data.other_mortality[i])

    def calculate_survival_second_age_group(self, uid_data: pd.DataFrame, age_id: int, prev_age_id: int):
        """The second youngest age group.
        """
         #find first and second age group positions
        i = np.where(uid_data.age_group_id==age_id)
        j = np.where(uid_data.age_group_id==prev_age_id)
        
        #check to make sure that MIR is within bounds
        if uid_data.mi_ratio[i] < max_mi_ratio_other_age_groups(uid_data.other_mortality[i],
                                                                uid_data.P_s[j]):
            #solve for P_c (probability of death due to cause)
            uid_data.P_c[i] = bisect(
                    partial(self.second_age_group_equation, 
                    other_mortality=uid_data.other_mortality[i],
                    mir=uid_data.mir[i],
                    P_s_1=uid_data.P_s[j]), 0.0, 1.0
                )
        else:
            #MIR is too large to properly solve (survival will go negative) assign survival of 0
            uid_data.P_c[i] = 1 - uid_data.other_mortality


        #solve for P_s (probability of abs survival in 1 year)
        uid_data.P_s[i] = 1 - (uid_data.P_c[i] + uid_data.other_mortality[i])

    def calculate_survival_other_age_group(self, 
                                            uid_data: pd.DataFrame, 
                                            age_id: int, 
                                            prev_age_id: int, 
                                            prev_prev_age_id: int):
        """Any other age group.
        """
         #find reference and two previous age group positions
        i = np.where(uid_data.age_group_id==age_id)
        j = np.where(uid_data.age_group_id==prev_age_id)
        k = np.where(uid_data.age_group_id==prev_prev_age_id)
        
        #check to make sure that MIR is within bounds
        if uid_data.mi_ratio[i] < max_mi_ratio_other_age_groups(uid_data.other_mortality[i],
                                                                uid_data.P_s[j],uid_data.P_s[k]):
            #solve for P_c (probability of death due to cause)
            uid_data.P_c[i] = bisect(
                    partial(self.other_age_group_equation, 
                    other_mortality=uid_data.other_mortality[i],
                    mir=uid_data.mir[i],
                    P_s_1=uid_data.P_s[j], P_s_2=uid_data.P_s[k]), 0.0, 1.0
                )

        else:
            #MIR is too large to properly solve (survival will go negative) assign survival of 0
            uid_data.P_c[i] = 1 - uid_data.other_mortality

        #solve for P_s (probability of abs survival in 1 year)
        uid_data.P_s[i] = 1 - (uid_data.P_c[i] + uid_data.other_mortality[i])

    ### CODING OF EQUATIONS ###
    @staticmethod
    def first_age_group_equation(other_mortality: float, mir: float, P_c: float):
        """Equation for the first age group relating MI ratio and survival

        Args: 
            P_c (float): probability of death in a single year due to cause
            other_mortality (float): probability of death due to other causes
            mir (float): mortality incidence ratio
            
        Returns:
                float: value of the difference between the two sides of the equation
        """
        P_s = 1 - (P_c + other_mortality)

        right_hand_side = P_c/(P_c + other_mortality)*(1-1/5*sum([P_s**i for i in range(1,5)]))

        return right_hand_side - mir

    @staticmethod
    def second_age_group_equation(other_mortality: float, 
                                mir: float,
                                P_c: float,
                                P_s_1: float):
        """Equation for the second age group relating MI ratio and survival

        Args: 
            P_c (float): probability of death in a single year due to cause
            other_mortality (float): probability of death due to other causes
            mir (float): mortality incidence ratio
            P_s_1 (float): probability of survival of the previous (youngest) age group
        
        Returns:
            float: value of the difference between the two sides of the equation
        """
        P_s = 1 - (P_c + other_mortality)

        right_hand_side = P_c/(P_c + other_mortality)*(1-1/5*sum([P_s**i for i in range(1,5)])+1/5*sum([P_s_1**i for i in range(1,5)])*(1-P_s**5))
        
        return right_hand_side - mir

    
    @staticmethod
    def other_age_group_equation(other_mortality: float, 
                                mir: float,
                                P_c: float,
                                P_s_1: float,
                                P_s_2: float):
        """Equation for the other age groups relating MI ratio and survival

        Args: 
            P_c (float): probability of death in a single year due to cause
            other_mortality (float): probability of death due to other causes
            mir (float): mortality incidence ratio
            P_s_1 (float): probability of survival of the previous age group
            P_s_2 (float): probability of survival of the two previous age group
        
        Returns:
            float: value of the difference between the two sides of the equation
        """
        P_s = 1 - (P_c + other_mortality)

        right_hand_side = P_c/(P_c + other_mortality)*(1-1/5*sum([P_s**i for i in range(1,5)])+1/5*sum([P_s_1**i for i in range(1,5)])*(1-P_s**5)+1/5*P_s_1**5*complicated_sigma(4,P_s,P_s_2))
        
        return right_hand_side - mir

    ### HELPERS ###
    def complicated_sigma(n: int, prob_survival: float, P_s_2: float):
    '''needs a better name but means to solve: \sum_{i=1}**n P_s_2**i(1-prob_survival**{5-i})'''
        value = 0
        for i in range(1,n):
            value = value + P_s_2**i*(1-prob_survival**(5-i))
        return value

    def max_mi_ratio_other_age_groups(other_mortality: float, P_s_1: float, P_s_2: float):
    '''The bound that MI ratio cannot exceed in this monotonically decreasing function of
     survival predicting MIR is P_s = 0. This is solved to be the following.'''
        max = (1 - other_mortality)*(1+1/5*sum([P_s_1**i for i in range(1,5)])+1/5*P_s_1**5*sum([P_s_2**i for i in range(1,4)]))
        return max

    def max_mi_ratio_second_age_group(other_mortality: float, P_s_1: float):
    '''The bound that MI ratio cannot exceed in this monotonically decreasing function of
     survival predicting MIR is P_s = 0. This is solved to be the following.'''
        max = (1 - other_mortality)*(1+1/5*sum([P_s_1**i for i in range(1,5)]))
        return max