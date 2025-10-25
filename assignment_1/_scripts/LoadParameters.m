% Load parameters
% NOTE: So far, I use every parameter. 
% Model parameters & settings

% Pre-defined parameters
pa              = 1;
pAlpha          = 0.3;
pTau            = 0.15;
pAvgz           = 1;
pA              = 1;
pr              = 0.04;
pBeta           = 0.96;

% Placeholder parameters
load('_results/ParametersPS.mat');
pSigmaz         = pSigmazPS;
pEta            = pEtaPS;
pb              = pbPS;
pChi            = pChiPS;

% Solving iterations parameters
pStepSize       = 0.9;
pMaxIter        = 2000;
pMaxWage        = 10;

% Target moments
mWeights        = eye(4);
vTargetMoments  = [0.33; 0.06; 0.25; 0.70];

% Parameter bounds
% 1. pChi
% 2. PEta
% 3. pb
% 4. pSigmaz
vUpperBound     = [Inf; Inf; Inf; Inf];
vLowerBound     = [0;0;0;0];