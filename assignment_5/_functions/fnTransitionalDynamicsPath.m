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

%% 3. Shooting algorithm
pT              = size(vAPath,1);
K0              = Results1.K;
KT              = Results2.K;
vKGuess         = ones(pT, 1)*K0;
fnShooting      = @(x) fnEulerResidualPath(x,vAPath,K0,KT, Parameters);
Solution        = fsolve(fnShooting,vKGuess);

%% 4. Save results
% Compute other vectors
vKp             = Solution;
vK              = [K0; vKp(1:end-1)];

% Initialise
vN              = zeros(size(vK));
vC              = vN;
vR              = vC;
vW              = vR;

% Other variables
for kkk = 1:1:size(vK,1)
    % Compute
    K           = vK(kkk);
    Kp          = vKp(kkk);
    [N,C,R,W]   = fnFindNTransitions(K,Kp,AT,Parameters);

    % Save the results
    vN(kkk)     = N;
    vC(kkk)     = C;
    vR(kkk)     = R;
    vW(kkk)     = W;
end 

% Final save
Results.vK      = vK;
Results.vN      = vN;
Results.vC      = vC;
Results.vR      = vR;
Results.vW      = vW;

end 

