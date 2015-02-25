MATLAB code for:
Born to Export: Understanding Export Growth in Bangladesh’s Apparel and Textiles Industry,
by Diaz de Astarloa, Eaton, Krishna, Roberts, Rodríguez-Clare, and Tybout.

Written by James Tybout and Bernardo Diaz de Astarloa @ PSU (2010-2015).


Summary
————————————————————————————————
The code solves a dynamic programming problem for firms’ buyers search intensity and simulates a sample of firms. It also allows for an indirect inference estimation routine of specified parameters.   


Details
————————————————————————————————
main_bte.m calls the main routines of the code in the following order:
  1. setparams.m: assigns values to main parameters and settings.
  2. solve_bte.m: finds an equilibrium of the search problem as a function of the parameters.
  3. trajec_bte.m: simulates the model for a sample of firms and computes moments that can be used in estimation.
  4. plot_*.m: routines to plot and generate figures of interest. 

calibration_bte.m and distance_bte.m allow for an indirect inference algorithm to estimate parameters of interest. In particular:
	- calibration_bte.m: calls a genetic algorithm routine to minimize a distance metric between data and model-based moments.
	- distance_bte.m: calls main_bte.m and computes the distance metric to be minimized.  


Brief description of subroutines
————————————————————————————————
~ solve_bte.m:
	- hazz_nospill.m: calculates posterior beliefs about product appeal. 
	- tauchen.m: computes transition densities and grids of exogenous variables.
	- tauchen_grid.m: constructs grids for exogenous variables. 
	- pi_tilda.m: computes expected payoff to a match.
	- policy_h.m: computes optimal search intensity at home (policy functions from FOC taking as given the output from pi_tilda payoff). 
	- policy_f_bte: analogous of policy_h.m for the foreign market. 

~ trajec_bte.m:
	- dis_traj.m: generates exogenous shocks trajectories and maps them to a grid.
	- brooks.m: computes Brooks tables.
	- quant_dur.m: duration of matched by quantiles of sales.
