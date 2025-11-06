% Loading parameters

%% 1. Standard parameters
Parameters.pZeta        = 0.01;
Parameters.pVarphi      = 0.80;
Parameters.pPhi         = 0.5;
Parameters.pEta         = 2;
Parameters.pb           = 0.2;
Parameters.pBarh        = 10;
Parameters.pr           = 0.04;
Parameters.pBeta        = 0.96;
Parameters.pSigmaz      = 0.05;
Parameters.pRho         = 0.90;
Parameters.pGammah      = 0.1;
Parameters.pGammaz      = 0.1;
Parameters.pGamma0      = 0.2;
Parameters.pT           = 50;

%% 2. Grids: Productivity, human capital, age, and wealth

% Tauchen productivity parameters
Parameters.pTauchenN    = 1.5;
Parameters.pTauchenM    = 11;

% Tauchen productiviy grid % INITIATE THE GRID INSIDE THE FUNCTION
[Grids.vGridZ,... 
    Grids.mTransitionZ] = fnTauchenLogNormal(Parameters);

% Age and human capital and age grids % INITIATE THE GRID INSIDE THE FUNCTION
Grids.vGridH            = linspace(0,Parameters.pBarh,Parameters.pBarh+1)';
Grids.vGridAge          = linspace(0,Parameters.pT,Parameters.pT+1)';

% Wealth grid parameters
Parameters.pMina        = 0;
Parameters.pMaxa        = 30;
Parameters.pNumGrida    = 60;

% Wealth grid  % INITIATE THE GRID INSIDE THE FUNCTION
Grids.vGrida            = fnWealthGridMMV(Parameters.pMina, Parameters.pMaxa, Parameters.pNumGrida);