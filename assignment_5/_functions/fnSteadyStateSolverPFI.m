function [Policy,EqVals]            = fnSteadyStateSolverPFI(A,Parameters,Grids)
    
%% 1. Unpack parameters & grids

% Parameters
pSigma              = Parameters.pSigma;
pAlpha              = Parameters.pAlpha;
pBeta               = Parameters.pBeta;
pChi                = Parameters.pChi;
pMu                 = Parameters.pMu;
pEta                = Parameters.pEta;
pDelta              = Parameters.pDelta;

% Grids
vGridK              = Grids.vWealthGrid;


%% 2. Iteration options
vK                  = vGridK;
vKNextGuess         = 1.05*vGridK;
vKNextNew           = zeros(size(vKNextGuess));
iError              = 10.0;
pErrorTol           = 1e-5;
pStepSize           = 0.9;
iIterNum            = 1;

%% 3. Define key functions

% For prices
Interest            = @(N,K)    pAlpha * K^(pAlpha-1) * N^(1-pAlpha) - pDelta;
Wage                = @(N,K)    (1-pAlpha) * K^(pAlpha) * N^(-pAlpha);
ConsImplied         = @(N,K)    (N^(1/pChi)*pEta / Wage(N,K))^(-1/pSigma);
Adjustment          = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
BCError             = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - ConsImplied(N,K) - Kp - Adjustment(K,Kp);
BCConsumption       = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - Kp - Adjustment(K,Kp);

% For Euler
Psi_1               = @(K,Kp) -pMu * (Kp-K) / K + Adjustment(K,Kp) / K;
Psi_2               = @(K,Kp) pMu * (Kp - K) / K;
ValueDerivative     = @(N,K,Kp) BCConsumption(N,K,Kp)^(-pSigma) * (1 + Interest(N,K) - Psi_1(K,Kp));
EulerError          = @(N,K,Kp,Np,Kpp) BCConsumption(N,K,Kp)^(-pSigma) * (1 + Psi_2(K,Kp)) - pBeta *  ValueDerivative(Np,Kp,Kpp);

tic;
%% 4. Start iterating 
while iError > pErrorTol
    for kkk = 1:1:size(vK,1)

        % A. Optimal Kp_new
        EulerErrorK         = @(x) fnEulerResidual(vK(kkk),x,vK,vKNextGuess,EulerError,BCError,ConsImplied,Interest,Wage)^2;
        Kp_new              = fminbnd(EulerErrorK,0,max(vK));

        % B. Update results 
        vKNextNew(kkk)      = Kp_new;
    end
    
    % C. Iteration items 
    vKNextGuess             = pStepSize*vKNextGuess + (1 - pStepSize)*vKNextNew;
    iError                  = max(abs(vKNextGuess-vKNextNew),[],"all");
    fprintf('Error:             %.3f \n', iError);
end 
Policy              = vKNextGuess;
toc;

%% 5. Find the steady state capital

% Start
K_ss                = 1;
iErrorK_ss          = 100;

while iErrorK_ss > 1e-10
    K_ss_next       = interp1(vK,Policy,K_ss,"pchip");
    iErrorK_ss          = abs(K_ss_next-K_ss);
    K_ss            = pStepSize*K_ss+ (1-pStepSize)*K_ss_next;
    fprintf('SUMMARY \n')
    fprintf('Error:             %.3f \n', iError);
    fprintf('Steady state K:    %.3f \n', K_ss);
end 

%% 6. Find the key equilibrium values

% What we know
EqVals.K                = K_ss;
[N_eq,C_eq,R_eq,W_eq]   = fnFindN(K_ss,K_ss,BCError,ConsImplied,Interest,Wage);
EqVals.N                = N_eq;
EqVals.C                = C_eq;
EqVals.R                = R_eq;
EqVals.W                = W_eq;
end
