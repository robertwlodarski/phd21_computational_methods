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

%% 2. Standard shock [A=1 to A=1.1]

Results         = fnAiyagari1994Transition(1.0,1.1,Parameters,Grids);
ResultsPE       = fnAiyagari1994TransitionPE(1.0,1.1,Parameters,Grids);

%% 3. A quick drop & immediate recovery [A=1 to A=0.95 then recovery]

Grids.vAPath    = Paths.v2;
Results2        = fnAiyagari1994Transition(1.0,1.0,Parameters,Grids);

%% 4. A quick drop & slow recovery [A=1 to A=0.95 then recovery]

Grids.vAPath    = Paths.v3;
Results3        = fnAiyagari1994Transition(1.0,1.0,Parameters,Grids);


