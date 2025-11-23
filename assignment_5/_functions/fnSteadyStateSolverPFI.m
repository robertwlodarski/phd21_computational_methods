function [vK]            = fnSteadyStateSolverPFI(A,Parameters,Grids)
    
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

%% 3. Define key functions
Interest            = @(N,K)    pAlpha * K^(pAlpha-1) * N^(1-pAlpha) - pDelta;
Wage                = @(N,K)    (1-pAlpha) * K^(pAlpha) * N^(-pAlpha);
ConsImplied         = @(N,K)    (N^(1/pChi)*pEta / Wage(N,K))^(-1/pSigma);
Adjustment          = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
BCError             = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - ConsImplied(N,K) - Kp - Adjustment(K,Kp);

%% 3. Start iterating 

while iError > pErrorTol
    for kkk = 1:1:size(vK,1)
        % A. Solve for today's N and C
        % Narrow down the function
        BCErrorK            = @(N) BCError(N,vK(kkk),vKNextGuess(kkk));
        LabourImp           = fzero(BCErrorK,1);
        ConsImp             = ConsImplied(LabourImp, vK(kkk));
        
        % B. Obtain tomorrow's K prime (Kpp)
        Kpp_bot                         = sum(vK<vKNextGuess(kkk));
        Kpp_bot(Kpp_bot<1)              = 1;
        Kpp_bot(Kpp_bot>=size(vK,1))    = size(vK,1)-1;
        Kpp_up                          = Kpp_bot+1;
        Bot_weight                      = (vKNextGuess(Kpp_up)-vKNextGuess(kkk)/(vKNextGuess(Kpp_up)-vKNextGuess(Kpp_bot)));
        Kpp                             = Bot_weight * vKNextGuess(Kpp_bot) + (1 - Bot_weight) * vKNextGuess(Kpp_up);

        % C. Obtain tomorrow's prices
        BCErrorK            = @(N) BCError(N,vKNextGuess(kkk),Kpp);
        Np                  = fzero(BCErrorK,1);
        Ip                  = Interest(Np,vKNextGuess(kkk));
        Wp                  = Wage(Np,vKNextGuess(kkk));
        
        % D. Optimal K_new 
        vKNextNew(kkk)      = (1 + Ip) * Kpp + Np * Wp - Adjustment(vKNextGuess(kkk),Kpp);
    end
    
    % E. Iteration items 
    vKNextGuess             = pStepSize*vKNextGuess + (1 - pStepSize)*vKNextNew;
    iError                  = max(abs(vKNextGuess-vKNextNew),[],"all");
end 


end