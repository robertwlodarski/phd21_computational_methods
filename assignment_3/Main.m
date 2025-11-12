% Aiyagari (1994)

%% 1. Prepare the environment
clear; 
close all;

% Add extra paths
addpath _functions/
addpath _scripts/

% Load parameters
LoadParameters;

%% 2. Simulate the model

fnSimulateAiyagari1994(Parameters,Grids)