%% Loading parameters

% Pre-set parameters
Parameters.pSigma       = 1.0;
Parameters.pRho          = 0.95;
Parameters.pSigmaEps    = 0.009;
Parameters.pAlpha       = 0.33;
Parameters.pBeta        = 0.99;
Parameters.pChi         = 1.0;
Parameters.pMu          = 0.60;
Parameters.pDelta       = 0.025;

% Temporary parameter
Parameters.pEta         = 1.0;


%% Grids 

% Grid options
Parameters.pGridAMin    = 1e-2;
Parameters.pGridAMax    = 25;
Parameters.pNumGridA    = 50;

% Compute grid
Grids.vWealthGrid       = linspace(Parameters.pGridAMin,Parameters.pGridAMax,Parameters.pNumGridA)';