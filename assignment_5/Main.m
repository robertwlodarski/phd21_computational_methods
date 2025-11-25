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

% Run the shock
ResultsJump             = fnTransitionalDynamicsJump(1,1.1,Paths.Path1,Parameters,Grids);
ResultsNegative         = fnTransitionalDynamicsPath(1,1,Paths.Path2,Parameters,Grids);