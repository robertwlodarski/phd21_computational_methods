%% 1. Loading parameters

% Pre-set parameters
Parameters.pRhoA        = 0.859;
Parameters.pSigmaA      = 0.014;
Parameters.pPhi         = 0.975;
Parameters.pEta         = 2.400;
Parameters.pAlpha       = 0.256;
Parameters.pGamma       = 0.640;
Parameters.pBeta        = 0.977;
Parameters.pDelta       = 0.069;

%% 2. TFP grid and Tauchen
Parameters.pTauchenN    = 3;
Parameters.pTauchenM    = 7;
% Grid
[Grids.vGridTFP,...
    Grids.mTransitionTFP]   = fnTauchenLogNormalA(Parameters);


