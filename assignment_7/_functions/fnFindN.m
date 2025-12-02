function [N,C,R,W]  = fnFindN(K,Kp,pA, Parameters)
    %% 1. Parameters & functions
    pSigma              = Parameters.pSigma;
    pAlpha              = Parameters.pAlpha;
    pChi                = Parameters.pChi;
    pMu                 = Parameters.pMu;
    pEta                = Parameters.pEta;
    pDelta              = Parameters.pDelta;

    % Functions
    Interest            = @(N,K)    pA * pAlpha * K^(pAlpha-1) * N^(1-pAlpha) - pDelta;
    Wage                = @(N,K)    pA * (1-pAlpha) * K^(pAlpha) * N^(-pAlpha);
    ConsImplied         = @(N,K)    (N^(1/pChi)*pEta / Wage(N,K))^(-1/pSigma);
    Adjustment          = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
    BCError             = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - ConsImplied(N,K) - Kp - Adjustment(K,Kp);

    % 1. Solve for N
    BCErrorKSq      = @(N) BCError(N,K,Kp)^2;
    N               = fminbnd(BCErrorKSq,1e-5,3-1e-5);

    % 2. Solve for other elements
    R               = Interest(N,K);
    W               = Wage(N,K);
    C               = (1+R)*K + N*W - Kp - Adjustment(K,Kp);
 
end 