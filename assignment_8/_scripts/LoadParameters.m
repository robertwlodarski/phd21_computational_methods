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

% Parameters
Parameters.pTauchenN    = 3;
Parameters.pTauchenM    = 7;

% Grid
[Grids.vGridTFP,...
    Grids.mTransitionTFP]   = fnTauchenLogNormalA(Parameters);

%% 3. Capital grid

% Grid settings
Parameters.pGridKMin    = 1e-5;
Parameters.pGridKMax    = 15;
Parameters.pNumGridK    = 50;

% Wealth grids
Grids.vGridK            = fnWealthGridMMV(Parameters.pGridKMin,Parameters.pGridKMax, Parameters.pNumGridK);

