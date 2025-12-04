%% 1. Loading parameters

% Pre-set parameters
Parameters.pSigma       = 1.0;
Parameters.pRhoA        = 0.950;
Parameters.pSigmaA      = 0.009;
Parameters.pRhoZ        = 0.90;
Parameters.pSigmaZ      = 0.05;
Parameters.pAlpha       = 0.33;
Parameters.pBeta        = 0.99;
Parameters.pChi         = 1.0;
Parameters.pDelta       = 0.025;
Parameters.pMu          = 0.60;
Parameters.pWealthMin   = 1;
Parameters.pA           = 1;

% Temporary parameter
Parameters.pEta         = 7.9;

%% 2. Labour endowment grid
Parameters.pTauchenN    = 3;
Parameters.pTauchenM    = 7;
% Grid
[Grids.vGridZ,...
    Grids.mTransitionZ]     = fnTauchenLogNormalZ(Parameters);
[Grids.vGridTFP,...
    Grids.mTransitionTFP]   = fnTauchenLogNormalA(Parameters);

%% 3. Other grids 
Parameters.pGridAMin    = Parameters.pWealthMin;
Parameters.pGridAMax    = 300;
Parameters.pNumGridA1   = 100;
Parameters.pNumGridA2   = 150;
% Wealth grids
Grids.vGridA1           = fnWealthGridMMV(Parameters.pGridAMin, Parameters.pGridAMax, Parameters.pNumGridA1);
Grids.vGridA2           = fnWealthGridMMV(Parameters.pGridAMin, Parameters.pGridAMax, Parameters.pNumGridA2);

