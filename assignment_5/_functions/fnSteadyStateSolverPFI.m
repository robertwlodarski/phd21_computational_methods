function [Results]            = fnSteadyStateSolverPFI(A,Parameters,Grids)
    
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
vKNextGuess         = vGridK;
vKNextNew           = zeros(size(vKNextGuess));
iError              = 10.0;
pErrorTol           = 1e-4;
pStepSize           = 0;

%% 3. Define key functions

% For prices
Interest            = @(N,K)    A * pAlpha * K^(pAlpha-1) * N^(1-pAlpha) - pDelta;
Wage                = @(N,K)    A * (1-pAlpha) * K^(pAlpha) * N^(-pAlpha);
ConsImplied         = @(N,K)    (N^(1/pChi)*pEta / Wage(N,K))^(-1/pSigma);
Adjustment          = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
BCError             = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - ConsImplied(N,K) - Kp - Adjustment(K,Kp);
BCConsumption       = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - Kp - Adjustment(K,Kp);

% For Euler
Psi_1               = @(K,Kp) -pMu * (Kp-K) / (K) + Adjustment(K,Kp) / K;
Psi_2               = @(K,Kp) pMu * (Kp - K) / K;
ValueDerivative     = @(N,K,Kp) BCConsumption(N,K,Kp)^(-pSigma) * (1 + Interest(N,K) - Psi_1(K,Kp));
EulerError          = @(N,K,Kp,Np,Kpp) BCConsumption(N,K,Kp)^(-pSigma) * (1 + Psi_2(K,Kp)) - pBeta *  ValueDerivative(Np,Kp,Kpp);

tic;
%% 4. Start iterating 
while iError > pErrorTol
    for kkk = 1:1:size(vK,1)

        % A. Optimal Kp_new
        opts                = optimoptions('fsolve','Display','none');
        EulerErrorK         = @(x) fnEulerResidual(vK(kkk),x,vK,vKNextGuess,Parameters,A);
        Kp_new              = fsolve(EulerErrorK,vK(kkk),opts);

        % B. Update results 
        vKNextNew(kkk)      = Kp_new;
    end
    
    % C. Iteration items 
    iError                  = max(abs(vKNextGuess-vKNextNew),[],"all");
    vKNextGuess             = pStepSize * vKNextGuess + (1 - pStepSize) * vKNextNew;
    fprintf('Error:             %.3f \n', iError);
end 
Policy              = vKNextGuess;
toc;

%% 5. Find the steady state capital

% Start
K_ss                = 1;
iErrorK_ss          = 100;

% Run the loop
while iErrorK_ss > 1e-10
    K_ss_next       = interp1(vK,Policy,K_ss,"pchip");
    iErrorK_ss      = abs(K_ss_next-K_ss);
    K_ss            = pStepSize*K_ss+ (1-pStepSize)*K_ss_next;
end 

%% 6. Find the key equilibrium values

% What we know
Results.K                   = K_ss;
[N_eq,C_eq,R_eq,W_eq]       = fnFindN(K_ss,K_ss,A,Parameters);
Results.N                   = N_eq;
Results.C                   = C_eq;
Results.R                   = R_eq;
Results.W                   = W_eq;
Results.Policy              = Policy;
end
