Tsunami Projekt: User Documentation
===================================

Installing and Running
----------------------

1. clone the repository with :code:`git clone ...` 
2. add and update the submodules with :code:`git submodule init` and :code:`git submodule update`
3. build with :code:`scons`
4. run the solver with :code:`./build/tsunami_lab NUMBER_OF_CELLS SOLVER_TYPE [-u "setup arg1 arg2"]` 
5. execute the tests with :code:`./build/tests`

The output of the Dam Break Problem is in :code:`/solutions`

bitte abändern, falsche codezeilen

Command line parameters when executing
--------------------------------------

:code:`NUMBER_OF_CELLS` = Number of cells into which the simulation is discretized (where number >= 1)
:code:`SOLVER_TYPE` = type of the solver which is used (:code:`ROE` or :code:`FWAVE`)
:code:`[-u "setup arg1 arg2"]` = choose between :code:`'DamBreak1d h_l h_r'`, :code:`'ShockShock1d h hu'` and :code:`'RareRare1d h hu'`, default is :code:`'DamBreak1d 10 5'`, args are to be input as floats
