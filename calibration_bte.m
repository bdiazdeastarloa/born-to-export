% -------------------------------------------------------------------------
% Born to export 
%
% calibration_bte: set up genetic algorithm estimation routine.
%
% Written by Bernardo Diaz de Astarloa @ PSU 2015.
% -------------------------------------------------------------------------

clear all
clc

fprintf('Born to Export.\n')
fprintf('\n')
fprintf('--------------------------------------------------------------\n')

tic;

codepath = pwd;
filename = [codepath,'/results/results_ga-estimation.mat'];
format long;

rng(82082);

% Initial condition (# of individuals x size of parameter vector).
% Parameters are: 
% 1. alpha (success, beta distribution)
% 2. beta (failure, beta distribution) 
% 3. 1/delta (Match separation rate)

pop = [0.716  3.161  0.33;...
       1.500  4.500  0.25;...
       3.000  6.000  0.15];

nvar = size(pop,2);
lbound = [0.1  2  0.1];          % lower bound for search
ubound = [4  9  1];              % upper bound for search

% Set GA algorithm options.
% Stop if no improvement for (secs)...
stall = 1*3600;
% Kill algorithm after (secs)...
kill = 20*3600;

% Allow for parallelization by setting 'UseParallel' to 'always'.
options = gaoptimset('Display','iter','PopulationSize',3,'Generations',100,...
'StallTimeLimit',stall,'TimeLimit',kill,'MutationFcn',@mutationadaptfeasible,...
'FitnessScalingFcn',@fitscalingrank,'InitialPopulation',pop,'UseParallel','always',...
'PlotFcns',@gaplotbestf,'EliteCount',0);

[Y,fval,exitflag,output,population,scores] = ga(@(Y) distance_bte(Y),nvar,[],[],[],[],lbound,ubound,[],options);  

ga_time = toc;

save filename;

