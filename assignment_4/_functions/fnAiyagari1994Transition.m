function [Results]          = fnAiyagari1994Transition(Start,End,Parameters,Grids)

%% 1. Unpacking
% Parameters
ParametersUsed          = Parameters;

% Grids
vGridZ                  = Grids.vGridZ;
mTransitionZ            = Grids.mTransitionZ;

%% 2. Compute the initial and terminal equilibria

% Terminal state
ParametersUsed.pA       = End;
ResultsT                = fnSolveAiyagari1994Iteration(ParametersUsed,Grids);

% Initial state
ParametersUsed.pA       = Start;
Results0                = fnSolveAiyagari1994Iteration(ParametersUsed,Grids);

%% 3. Set-up
% Iterations business: Convergence, acceleration, and GE
iWeightOld              = 0.9;
iErrorGE                = 10;
iTolGE                  = 1e-8;
iTolDist                = 1e-8;
iAccelerationInterv     = 20;
iAccelerationStart      = 30;
iIterNumGE              = 1;

%% 4. Stationary dist. of labour allocation
iErrorZDist             = 10;
iTolZDist               = 1e-8;
vZDist                  = ones(size(vGridZ,1),1) / sum( ones(size(vGridZ,1),1));
while iErrorZDist>iTolZDist
    vZDistNext          = mTransitionZ' * vZDist;
    iErrorZDist         = abs(vZDistNext-vZDist);
    vZDist              = vZDistNext; 
end
vLabSupply              = vGridZ'*vZDist;

%% 5. TRANSITIONAL DYNAMICS

end 