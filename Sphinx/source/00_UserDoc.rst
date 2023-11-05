0. Tsunami Projekt: User Documentation
======================================

Installing and Running
----------------------

1. clone the repository with :code:`git clone https://github.com/darafsa/tsunami_lab.git` 
2. add and update the submodules with :code:`git submodule init` and :code:`git submodule update` 
3. build with :code:`scons` 
4. run the solver with :code:`./build/tsunami_lab CELLS SOLVER SETUP [height] [velocity]` 
5. execute the tests with :code:`./build/tests` 

Command line parameters when executing
--------------------------------------

| :code:`CELLS` = Number of cells in x-direction (where number >= 1) 
| :code:`SOLVER` = Type of the solver which (:code:`ROE` or :code:`FWAVE`) 
| :code:`SETUP` = Setup to use (:code:`DAMBREAK`, :code:`RARE` or :code:`SHOCK`) 
| :code:`[height, velocity]` (optional) = The height and velocity to use for RareRare and ShockShock Setup 
