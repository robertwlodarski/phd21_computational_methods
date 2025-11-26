%% 1. Loading parameters

% Pre-set parameters
Parameters.pSigma       = 1.0;
Parameters.pRho         = 0.95;
Parameters.pSigmaEps    = 0.009;
Parameters.pAlpha       = 0.33;
Parameters.pBeta        = 0.99;
Parameters.pChi         = 1.0;
Parameters.pMu          = 0.60;
Parameters.pDelta       = 0.025;

% Temporary parameter
% Parameters.pEta         = 1.0;


%% 2. Grids 

% Grid options
Parameters.pGridAMin    = 1e-2;
Parameters.pGridAMax    = 25;
Parameters.pNumGridA    = 50;

% Compute grid
Grids.vWealthGrid       = linspace(Parameters.pGridAMin,Parameters.pGridAMax,Parameters.pNumGridA)';

% Tauchen grid for A
Parameters.pTauchenN    = 3; 
Parameters.pTauchenM    = 7;
[Grids.vGridA,...
    Grids.mTransitionA] = fnTauchenLogNormal(Parameters);

%% 3. Paths for transition

% Set initial path
Paths.pT                = 200;
Paths.Path1             = ones(Paths.pT,1)*1.1;

% Set up the logarithmic path
Paths.Path2             = zeros(Paths.pT,1);
% Temporary A path
Paths.Path2(1)          = 0.95;
for ttt=2:1:Paths.pT
    Paths.Path2(ttt)=Paths.Path2(ttt-1)^(0.95);
end
clear ttt;
