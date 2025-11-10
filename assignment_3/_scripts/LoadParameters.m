% Load parameters 

%% Pre-set parameters
Parameters.pRho         = 0.9;
Parameters.pSigma       = 0.2;
Parameters.pAlpha       = 0.36;
Parameters.pBeta        = 0.96;
Parameters.pDelta       = 0.08;
Parameters.pA           = 1;

%% Labour endowment grid
Parameters.pTauchenN    = 3;
Parameters.pTauchenM    = 7;
% Grid
[Grids.vGridZ,...
    Grids.mTransitionZ] = fnTauchenLogNormal(Parameters);