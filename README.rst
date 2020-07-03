=================
Survival Analysis
=================


This repository is created for the applications from cancer team.
The goal of the program is to infer the survival rate simply from
mortality incidence ratio of the diease and the background mortality.


Install
-------
The package require ``python=3.7`` and preliminary packages, Numpy, Scipy,
and Pytest. If want to run the examples, Pandas, Jupyter and Maplotlib are also required.

To install, simply clone the repository and install,

.. code-block:: shell

   git clone https://github.com/ihmeuw-msca/Survival.git && cd Survival
   python setup.py install


Example
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

For the result, ``survival_rate`` is a dictionary with ``survival_rate['abs']``
store the absolute survival rate and ``survival_rate['rel']`` store the relative
survival rate. And the argument ``num_years=5`` indicates that we compute the
5 years survival rate.

For more detailed examples please check `examples <https://github.com/ihmeuw-msca/Survival/blob/master/examples>`_.