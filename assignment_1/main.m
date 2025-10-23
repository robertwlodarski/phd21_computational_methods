% PhD21 Computational methods: Assignment 1
% Rob Wlodarski
% October 2025

%% 1. Prepare the environment
clear; 
close all;

% Add extra paths
addpath _functions/
addpath _scripts/

%% 2. Solve the model

% Load parameters
parameters;

% Calibrate the model (standard)
calibrate;

% Calibrate the model (alternative options)
calibrate_pattern_search;