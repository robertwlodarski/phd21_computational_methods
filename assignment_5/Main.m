% PhD21 Computational methods: Assignment 5
% Rob Włodarski
% November 2025

clear;
clc;

% Add key paths
addpath _functions/
addpath _scripts/

%% 1. Prepare parameters
LoadParameters;

%% 2. Transitional dynamics

% A. Dry runs for A=1 and A=1.1.
Results1        = fnSteadyStateSolverPFI(1,Parameters,Grids);
Results2        = fnSteadyStateSolverPFI(1.1,Parameters,Grids);