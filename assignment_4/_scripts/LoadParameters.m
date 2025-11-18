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

%% 4. Transitional dynamics related
Parameters.pT           = 200;

% Path 1 
Paths.v1                = repmat(1.1,Parameters.pT,1);
Paths.v1(1)             = 1;
Grids.vAPath            = Paths.v1;

% Path 2
Paths.v2                = ones(Parameters.pT,1);
Paths.v2(1)             = 0.95;

% Path 3
Paths.pRho              = 0.90;
Paths.v3                = zeros(Parameters.pT,1);
Paths.v3(1)             = 0.95;
for ttt = 2:1:Parameters.pT
    Paths.v3(ttt)       = Paths.v3(ttt-1)^(Paths.pRho);
end
