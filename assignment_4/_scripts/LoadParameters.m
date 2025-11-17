% Load parameters 

%% 1. Pre-set parameters
Parameters.pRho         = 0.9;
Parameters.pSigma       = 0.2;
Parameters.pAlpha       = 0.36;
Parameters.pBeta        = 0.96;
Parameters.pDelta       = 0.08;
Parameters.pA           = 1;

%% 2. Labour endowment grid
Parameters.pTauchenN    = 3;
Parameters.pTauchenM    = 7;
% Grid
[Grids.vGridZ,...
    Grids.mTransitionZ] = fnTauchenLogNormal(Parameters);

%% 3. Other grids 
Parameters.pGridAMin    = 0;
Parameters.pGridAMax    = 150;
Parameters.pNumGridA1   = 50;
Parameters.pNumGridA2   = 100;
% Wealth grids
Grids.vGridA1           = fnWealthGridMMV(Parameters.pGridAMin, Parameters.pGridAMax, Parameters.pNumGridA1);
Grids.vGridA2           = fnWealthGridMMV(Parameters.pGridAMin, Parameters.pGridAMax, Parameters.pNumGridA2);

%% 4. Simulation related
Parameters.pN           = 10000;
Parameters.pT           = 50;