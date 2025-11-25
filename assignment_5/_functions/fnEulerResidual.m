function [Res]           = fnEulerResidual(K,Kp,vK,vKNextGuess,Parameters,pA)

%% 1. Parameters & functions
pSigma              = Parameters.pSigma;
pAlpha              = Parameters.pAlpha;
pMu                 = Parameters.pMu;
pDelta              = Parameters.pDelta;
pBeta               = Parameters.pBeta;

% For prices
Interest            = @(N,K)    pA * pAlpha * K^(pAlpha-1) * N^(1-pAlpha) - pDelta;
Wage                = @(N,K)    pA * (1-pAlpha) * K^(pAlpha) * N^(-pAlpha);
Adjustment          = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
BCConsumption       = @(N,K,Kp) (1 + Interest(N,K))*K + Wage(N,K)*N - Kp - Adjustment(K,Kp);

% For Euler
Psi_1               = @(K,Kp) -pMu * (Kp-K) / (K) - Adjustment(K,Kp) / K;
Psi_2               = @(K,Kp) pMu * (Kp - K) / K;
ValueDerivative     = @(N,K,Kp) BCConsumption(N,K,Kp)^(-pSigma) * (1 + Interest(N,K) - Psi_1(K,Kp));
EulerError          = @(N,K,Kp,Np,Kpp) BCConsumption(N,K,Kp)^(-pSigma) * (1 + Psi_2(K,Kp)) - pBeta *  ValueDerivative(Np,Kp,Kpp);

%% 2. Find contemporaneous N, W, C and R
[N,C,~,~]               = fnFindN(K,Kp,pA,Parameters);

% Break if negative consumption
if C <= 1e-6
    % Return a massive error that scales with how badly we violated the constraint.
    % This guides fminbnd back toward feasible K'.
    Res = 1e10 + (1e-6 - C)^2;
    return;
end

%$ 3. Fast interpolation of K''
Kpp_bot                        = sum(vK<=Kp);
Kpp_bot(Kpp_bot<1)             = 1;
Kpp_bot(Kpp_bot>=size(vK,1))   = size(vK,1)-1;
Kpp_up                         = Kpp_bot+1;
Bot_weight                     = (vK(Kpp_up)-Kp)/(vK(Kpp_up)-vK(Kpp_bot));
Kpp                            = Bot_weight * vKNextGuess(Kpp_bot) + (1 - Bot_weight) * vKNextGuess(Kpp_up);

%% 4. Get tomorrow's elements
[Np,~,~,~]              = fnFindN(Kp,Kpp,pA,Parameters);

%% 5. Get the residual
Res                     = EulerError(N,K,Kp,Np,Kpp);

end 