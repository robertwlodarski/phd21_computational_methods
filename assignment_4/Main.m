% Assignment 4: Aiyagari (1994) w/ transitional dynamics

%% 1. Prepare the environment

% Clearing
clear; 
close all;
clear global;

% Add paths
addpath _functions\
addpath _scripts\

% Load baseline parameters
LoadParameters;

%% 2. Questions (a) & (b) 

Results         = fnAiyagari1994Transition(1.0,1.1,Parameters,Grids);

