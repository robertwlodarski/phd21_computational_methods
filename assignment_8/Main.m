% PhD21 Computational methods: Assignment 8
% Rob Włodarski
% January 2026

clear;
clc;

% Add key paths
addpath _functions/
addpath _scripts/

%% 1. Prepare parameters
LoadParameters;

%% 2. Solve the steady state model
ResultSS            = fnSteadyState(Parameters);

%% 3. Run the uncertainty equilibrium
Results             = fnUncertaintyEquilibrium(Parameters,Grids, ResultSS);
