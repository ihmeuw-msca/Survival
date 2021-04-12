##################################################
# Author: Franny Dean
# Purpose: Age group extension to survival model v2
# Date: 9/21/2020
#
##################################################

from typing import Dict
from numbers import Number
from functools import partial
import pandas as pd
import numpy as np 
from scipy.optimize import bisect

class AgeSurvivalModel2:
    """ Mortality incidence ratio used to approximate age specific survival rate
    """

    def __init__(self, inputs: pd.DataFrame):
        """Constructor of the model.

            Arguments:
                inputs: data frame that has the required columns
        """
        required = ['age_group_id','location_id','other_mortality','sex_id','acause','year_id','deaths','cases']

        self.inputs = pd.DataFrame(inputs[required])

        assert [column for column in required if(column in self.inputs.columns.tolist())]

        #no NAs
        self.inputs = self.inputs.dropna()

        #length
        self.num_points = len(self.inputs)
    
        #check inputs variable values
        self._check_inputs()

        #place holders
        self.inputs['P_c'] = np.full(self.num_points, np.nan)
        self.inputs['P_s'] = np.full(self.num_points, np.nan)
        self.inputs['abs_survival'] = np.full(self.num_points, np.nan)
        self.inputs['rel_survival'] = np.full(self.num_points, np.nan)

        #create UID column on location-year-sex-cause
        self.inputs['UID'] = self.inputs.apply(
                            lambda row: f"{row['location_id']} {row['year_id']} {row['sex_id']} {row['acause']}", axis=1
                            )
        #self.inputs['UID'] = str(self.inputs['location_id'])+str(self.inputs['year_id'])+str(self.inputs['sex_id'])+str(self.inputs['acause'])

        #create df dictionaries
        self.input_dfs = dict()
        
        for uid in self.inputs['UID'].unique():
            uid_data = self.inputs[self.inputs['UID']==uid].copy()
            #sort by age group (ascending)
            uid_data.sort_values(by=['age_group_id'], inplace=True)
            #import pdb; pdb.set_trace()
            uid_data.reset_index(inplace=True, drop=True)
            self.input_dfs[uid] = uid_data

    def combine_outputs(self):
        '''combines the UID data frames'''
        return pd.concat(list(self.input_dfs.values()),axis=0).reset_index(drop=True)

    def _check_inputs(self):
        """Check the values of the inputs for the class for suitability.
        """

        #want other_mortality less than 1, greater than zero
        assert (self.inputs['other_mortality']>=0).all()
        assert (self.inputs['other_mortality']<1).all()

    def compute_P_s(self):
        """Main function which computes the probablity of survival P_s for each
            age group in sequence.
        """
        for uid, data in self.input_dfs.items():
            #compute base case age groups
            self.calculate_survival_first_age_group(data)
            self.calculate_survival_second_age_group(data)

            #compute remaining age groups
            #import pdb; pdb.set_trace()
            for i in range(2,len(data.age_group_id)):
                self.calculate_survival_other_age_group(data, i, i-1, i-2)

    def compute_rel_survival(self):
        """Function for getting the 5-year relative survival.
        """
        for uid, data in self.input_dfs.items():
            #catch to make sure P_s has been calculated
            if np.isnan(data['P_s']).all():
                print("Computing survival first.")
                self.compute_P_s()

            #calculate survival directly from P_s
            for i in range(0,len(data)):
                data.loc[i,'abs_survival'] = data.P_s.tolist()[i]
                data.loc[i,'rel_survival'] = data.abs_survival.tolist()[i]/(1-data.other_mortality.tolist()[i]**5)

    ### SOLVE EQUATIONS ###
    def calculate_survival_first_age_group(self, uid_data: pd.DataFrame):
        """The youngest age group (below which there is assumed to be no incidence).
        """
        #first age group
        i=0

        #check to make sure that MIR is within bounds, needs to be under 1 for first age group
        if uid_data.cases.tolist()[i]==0: #no cases/infinite MIR, set probability of death by cancer at 0
            uid_data.loc[i,'P_c'] = 0

        else:
            if uid_data.deaths.tolist()[i]/uid_data.cases.tolist()[i] < 1:
            
                #import pdb; pdb.set_trace()
                #solve for P_c (probability of death due to cause)
                uid_data.loc[i,'P_c'] = bisect(
                    partial(self.first_age_group_equation, 
                    other_mortality=uid_data.other_mortality.tolist()[i],
                    deaths=uid_data.deaths.tolist()[i], cases=uid_data.cases.tolist()[i]), -100.0, 100.0)
            else:
                #MIR is too large to properly solve (survival will go negative) assign survival of 0
                uid_data.loc[i,'P_c']  = 1 - uid_data.other_mortality.tolist()[i]


        #solve for P_s (probability of abs survival in 1 year)
        uid_data.loc[i,'P_s']  = 1 - (uid_data.P_c.tolist()[i] + uid_data.other_mortality.tolist()[i])

    def calculate_survival_second_age_group(self, uid_data: pd.DataFrame):
        """The second youngest age group.
        """
        #first and second age group positions
        i = 1
        j = 0
        
        #check to make sure that MIR is within bounds
        if uid_data.cases.tolist()[i]==0: #no cases/infinite MIR, set probability of death by cancer at 0
            uid_data.loc[i,'P_c'] = 0

        else:
            the_max = 2 #self.max_mi_ratio_second_age_group(uid_data.other_mortality.tolist()[i], uid_data.P_s.tolist()[j])
            if uid_data.deaths.tolist()[i]/uid_data.cases.tolist()[i] <+ the_max:
                #solve for P_c (probability of death due to cause)
                uid_data.loc[i,'P_c']  = bisect(
                    partial(self.second_age_group_equation, 
                    other_mortality=uid_data.other_mortality.tolist()[i],
                    deaths=uid_data.deaths.tolist()[i], cases=uid_data.cases.tolist()[i],
                    cases_prev=uid_data.cases.tolist()[j],
                    P_s_1=uid_data.P_s.tolist()[j]), -100.0, 100.0)
            else:
                #MIR is too large to properly solve (survival will go negative) assign survival of 0
                uid_data.loc[i,'P_c']  = 1 - uid_data.other_mortality.tolist()[i]


        #solve for P_s (probability of abs survival in 1 year)
        uid_data.loc[i,'P_s']  = 1 - (uid_data.P_c.tolist()[i] + uid_data.other_mortality.tolist()[i])

    def calculate_survival_other_age_group(self, 
                                            uid_data: pd.DataFrame, 
                                            i: int, 
                                            j: int, 
                                            k: int):
        """Any other age group.
        """
        
        #check to make sure that MIR is within bounds
        if uid_data.cases.tolist()[i]==0: #no cases/infinite MIR, set probability of death by cancer at 0
            uid_data.loc[i,'P_c'] = 0

        else:
            the_max = 2 #self.max_mi_ratio_other_age_groups(uid_data.other_mortality.tolist()[i], uid_data.P_s.tolist()[j],uid_data.P_s.tolist()[k])
            if uid_data.deaths.tolist()[i]/uid_data.cases.tolist()[i]  <= the_max:
                #solve for P_c (probability of death due to cause)
                uid_data.loc[i,'P_c']  = bisect(
                    partial(self.other_age_group_equation, 
                    other_mortality=uid_data.other_mortality.tolist()[i],
                    deaths=uid_data.deaths.tolist()[i], cases=uid_data.cases.tolist()[i],
                    cases_prev=uid_data.cases.tolist()[j], cases_2_prev=uid_data.cases.tolist()[k],
                    P_s_1=uid_data.P_s.tolist()[j], P_s_2=uid_data.P_s.tolist()[k]), -100.0, 100.0)

            else:
                #MIR is too large to properly solve (survival will go negative) assign survival of 0
                uid_data.loc[i,'P_c']  = 1 - uid_data.other_mortality.tolist()[i]

        #solve for P_s (probability of abs survival in 1 year)
        uid_data.loc[i,'P_s'] = 1 - (uid_data.P_c.tolist()[i] + uid_data.other_mortality.tolist()[i])

    ### CODING OF EQUATIONS ###
    @staticmethod
    def first_age_group_equation(P_c: float, other_mortality: float, deaths: float, cases: float):
        """Equation for the first age group relating MI ratio and survival

        Args: 
            P_c (float): probability of death in a single year due to cause
            other_mortality (float): probability of death due to other causes
            mir (float): mortality incidence ratio
            
        Returns:
                float: value of the difference between the two sides of the equation
        """

        return cases*(P_c+other_mortality) - deaths

    @staticmethod
    def second_age_group_equation(P_c: float, other_mortality: float, 
                                deaths: float, cases: float,
                                cases_prev: float,
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

        right_hand_side = cases*(P_c+other_mortality) + cases_prev*P_s_1*(P_c+other_mortality)

        return right_hand_side - deaths

    @staticmethod
    def other_age_group_equation(P_c: float, other_mortality: float, 
                                deaths: float, cases: float,
                                cases_prev: float, cases_2_prev: float,
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

        right_hand_side = cases*(P_c+other_mortality) + cases_prev*P_s_1*(P_c+other_mortality) + cases_2_prev*P_s_1*P_s_2*(P_c+other_mortality)
        return right_hand_side - deaths

    ### HELPERS ###
    #@staticmethod
    #def max_mi_ratio_other_age_groups(other_mortality: float, P_s_1: float, P_s_2: float):
    #   '''The bound that MI ratio cannot exceed in this monotonically decreasing function of
    #    survival predicting MIR is P_s = 0. This is solved to be the following.'''
            
    #    the_max = (1 - other_mortality)*(1+1/5*sum([P_s_1**i for i in range(1,5)])+1/5*P_s_1**5*sum([P_s_2**i for i in range(1,4)]))
    #   return the_max

    #@staticmethod
    #def max_mi_ratio_second_age_group(other_mortality: float, P_s_1: float):
    #    '''The bound that MI ratio cannot exceed in this monotonically decreasing function of
    #    survival predicting MIR is P_s = 0. This is solved to be the following.'''
            
    #    the_max = (1 - other_mortality)*(1+1/5*sum([P_s_1**i for i in range(1,5)]))
    #    return the_max