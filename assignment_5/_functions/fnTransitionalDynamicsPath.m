function [Results]          = fnTransitionalDynamicsPath(A0,AT,vAPath,Parameters,Grids)

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
if AT == A0
    Results2    = Results1;
else
    Results2        = fnSteadyStateSolverPFI(AT,Parameters,Grids);
end

% Temporary A path
vAPath          = zeros(100,1);
vAPath(1)       = 0.95;
for ttt=2:1:100
    vAPath(ttt)=vAPath(ttt-1)^(0.95);
end


%% 3. Shooting algorithm
pT              = size(vAPath,1);
K0              = Results1.K;
KT              = Results2.K;
vKGuess         = ones(pT, 1)*K0;
fnShooting      = @(x) fnEulerResidualPath(x,vAPath,K0,KT, Parameters);
Solution        = fsolve(fnShooting,vKGuess);


end 