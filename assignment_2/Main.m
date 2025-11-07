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
Simulations     = fnSimulationsSolver(Parameters,Grids,20000);
toc;

% Plotting