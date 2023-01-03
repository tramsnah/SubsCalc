Python scripts to do Geertsma/van Opstal based calculations
===========================================================
Contains shapefile-based subsidence bowl calculation.
Additional utilities are provided to compute time dependence.
This includes a rudimentary pvt module to convert production history
to pressure history.
The general assumption is:
	subs(x,y,t) = b(x,y)*f(t)
for an individual bowl. Multiple bowls, with different time
dependencies can be added.

Folder			Purpose
---------------------------------------------------------------------------------------
py_opstal		Python script library
working			'Live' folder (contains lots of garbage as well as useful stuff)
subs_XXX		Computations for field XXX

Relative folder locations are used in the scripts so should not be changed.

To test the code, run:
py_opstal/vanopstal.py
py_opstal/subscalc.py
py_opstal/profile_utils.py
py_opstal/pvtcorrelation.py


