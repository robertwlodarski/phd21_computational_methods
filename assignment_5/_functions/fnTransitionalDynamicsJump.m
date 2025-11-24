function [Results]          = fnTransitionalDynamicsJump(A0,A1,Parameters,Grids)

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
% Useful functions
A               = A1;
Interest        = @(N,K)    A * pAlpha * K^(pAlpha-1) * N^(1-pAlpha) - pDelta;
Wage            = @(N,K)    A * (1-pAlpha) * K^(pAlpha) * N^(-pAlpha);
ConsImplied     = @(N,K)    (N^(1/pChi)*pEta / Wage(N,K))^(-1/pSigma);
Adjustment      = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
BCError         = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - ConsImplied(N,K) - Kp - Adjustment(K,Kp);

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

% Run the loop
for ttt = 1:1:pT
    
    % A. Carry previous period's result
    if ttt>1
        vK(ttt) = vKNext(ttt-1); 
    else 
    end
    
    % B. Find Kp
    vKNext(ttt) = interp1(vGridK, vPolicy, vK(ttt), 'pchip');

    % C. Find N et al.
    [N,C,R,W]   = fnFindN(vK(ttt),vKNext(ttt),BCError,ConsImplied,Interest,Wage);
    
    % D. Save results 
    vN(ttt)     = N;
    vC(ttt)     = C;
    vR(ttt)     = R;
    vW(ttt)     = W;

    % E. Break if no meaningul change
    % Set change
    Distance      = abs(vK(ttt) - Results2.K);
    % Break the loop if not needed
    if Distance < pTolDistance
        break
    else
        continue
    end 
end 

% Save observations without NaNs
Results.vTime   = sum(~isnan(vK));
Results.vK      = vK(~isnan(vK));
Results.vN      = vN(~isnan(vK));
Results.vC      = vC(~isnan(vK));
Results.vR      = vR(~isnan(vK));
Results.vW      = vW(~isnan(vK));

end 