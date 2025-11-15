% Aiyagari (1994)

%% 1. Prepare the environment
clear; 
close all;

% Add extra paths
addpath _functions/
addpath _scripts/

% Load parameters
LoadParameters;

%% 2. Solve the model using the iteration method
Results             = fnSolveAiyagari1994Iteration(Parameters,Grids);

%% 3. Simulations
Simulations         = fnGenerateSimulations(1000,10000,Results, Parameters, Grids);