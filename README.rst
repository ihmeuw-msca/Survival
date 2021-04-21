=================
Survival Analysis
=================


This repository is created for the applications from cancer team. The goal of the program is to infer the survival rate simply from mortality incidence ratio of the diease and the background mortality.


Install
-------
The package require ``python=3.7`` and preliminary packages, Numpy, Scipy,
and Pytest. If want to run the examples, Pandas, Jupyter and Maplotlib are also required.

To install, simply clone the repository and install,

.. code-block:: shell

   git clone https://github.com/ihmeuw-msca/Survival.git && cd Survival
   python setup.py install


Base Model Example
-------
After installation, to use the code, we need data for ``mir`` (mortality incidence ratio),
``other_mortality`` and optional ``disease_period``, for more details, please
check the `docstring <https://github.com/ihmeuw-msca/Survival/blob/master/src/survival/model.py#L20-L28>`_.
And we could create object and do the computation.

.. code-block:: python

   from survival import MIRModel

   model = MIRModel(mir, other_mortality, disease_period=5)
   model.compute_excess_mortality()
   
   survival_rate = model.get_survival_rate(num_years=5)

For the result, ``survival_rate`` is a dictionary with ``survival_rate['abs']`` store the absolute survival rate and ``survival_rate['rel']`` store the relative survival rate. And the argument ``num_years=5`` indicates that we compute the 5 years survival rate.

For more detailed examples please check `examples <https://github.com/ihmeuw-msca/Survival/blob/master/examples>`_.

Age Extension Model Example
-------
The GBD study is interested in predicting cancer burden in five year age bins. The age group extension model incorporates the fact that individuals who are diagnosed in one age group may die in another. In the model, the first age group uses the base MIR model to obtain an estimate of survival. The next oldest age groups incorporates survivors from the youngest age group. Older age groups use include the survival of the two previous age groups. A period of 10 years of excess mortality is assumed. The derivation of the model has been written up `_here <https://www.overleaf.com/read/hxfwhvsmmtnb>_`_. 

To use the code, we need a data frame with the following columns as listed `_here<https://github.com/ihmeuw-msca/Survival/blob/48a072a19544babfa204c443fccaa37d2babbc77/src/survival/age_group_extension.py#L31>_`_:

``[age_group_id','location_id','other_mortality','sex_id','acause','year_id','mi_ratio']``

.. code-block:: python
   from survival import AgeSurvivalModel
  
   model = AgeSurvivalModel(df_inputs)
   model.compute_P_s()

   model.compute_n_year_survival(num_years=5)
   
   df_outputs = model.combine_outputs()
   
The result ``df_outputs`` is the 5 year relative and absolute survival estimates from the model appended to the input data frame.

For more detailed jupyter notebook example uses please check `examples <https://github.com/ihmeuw-msca/Survival/blob/48a072a19544babfa204c443fccaa37d2babbc77/Running%20Age%20Group%20Extension.ipynb>`_.
