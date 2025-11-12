% PhD21 Computational methods: Assignment 2
% Rob Wlodarski
% October 2025

%% 1. Prepare the environement
clear;
close all;

% Add paths
addpath '_functions'; 
addpath '_scripts'; 
addpath '_plotting';

% Load parameters
LoadParameters;

% Simulate the model
pSimulationsNumber  = 10000;
tic;
[~,~,VFs.mW,VFs.mN] = fnValueFunctionMatrices(Parameters,Grids);
toc;
Simulations         = fnSimulationsSolver(Parameters,Grids,pSimulationsNumber,VFs);

% Plotting
StandardPlots;