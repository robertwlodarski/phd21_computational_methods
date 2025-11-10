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
tic;
pSimulationsNumber  = 20000;
[~,~,VFs.mW,VFs.mN] = fnValueFunctionMatrices(Parameters,Grids);
Simulations         = fnSimulationsSolver(Parameters,Grids,pSimulationsNumber,VFs);
toc;

% Plotting
StandardPlots;