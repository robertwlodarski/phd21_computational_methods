function [Results]          = fnTransitionalDynamicsPath(A0,A1,vAPath,Parameters,Grids)

%% 1. Parameters & grids
pSigma              = Parameters.pSigma;
pAlpha              = Parameters.pAlpha;
pChi                = Parameters.pChi;
pMu                 = Parameters.pMu;
pEta                = Parameters.pEta;
pDelta              = Parameters.pDelta;

% Grids
vGridK              = Grids.vWealthGrid;

%% 2. Run the steady state results
Results1        = fnSteadyStateSolverPFI(A0,Parameters,Grids);
Results2        = fnSteadyStateSolverPFI(A1,Parameters,Grids);

%% 3. Iterations
% Set time & tolerance for change
pT              = size(vAPath,1);
pTolDistance    = 1e-5;

% Prepare key items
vPolicy         = Results2.Policy;
vKNext          = nan(pT,1);
vC              = vKNext;
vW              = vKNext;
vR              = vKNext;
vN              = vKNext;
vK              = vKNext;
vK(1)           = Results1.K;
VKGuess         = linspace(Results1.K,Results2.K,pT);



% Save observations without NaNs
Results.vTime   = sum(~isnan(vK));
Results.vK      = vK(~isnan(vK));
Results.vN      = vN(~isnan(vK));
Results.vC      = vC(~isnan(vK));
Results.vR      = vR(~isnan(vK));
Results.vW      = vW(~isnan(vK));

end 