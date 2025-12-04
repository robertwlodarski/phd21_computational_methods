% PhD21 Computational methods: Assignment 7
% Rob Włodarski
% December 2025

clear;
clc;

% Add key paths
addpath _functions/
addpath _scripts/

%% 1. Prepare parameters
LoadParameters;

%% 2. Solve the model
Result0             = fnSolveSSIteration(Parameters,Grids);