function [vResiduals]        = fnEulerResidualPath(vKGuess,vAPath,K0,KT, Parameters)

%% 1. Parameters
% Unpack parameters
pSigma              = Parameters.pSigma;
pBeta               = Parameters.pBeta;
pMu                 = Parameters.pMu;

% For Euler
Adjustment          = @(K,Kp)   pMu / 2 * ((Kp - K)/K)^2*K;
Psi_1               = @(K,Kp) -pMu * (Kp-K) / (K) - Adjustment(K,Kp) / K;
Psi_2               = @(K,Kp) pMu * (Kp - K) / K;

% Other key items
pT              = size(vKGuess,1);
vResiduals      = zeros(size(vKGuess));
vKNext          = [K0;vKGuess; KT];
vAPathExpanded  = [vAPath; vAPath(end)];

%% 2. Run the loop
for ttt = 1:1:pT
    
    % A. K, Kp, Kpp, A, and Ap
    K                   = vKNext(ttt);
    Kp                  = vKNext(ttt+1);
    Kpp                 = vKNext(ttt+2);
    A                   = vAPathExpanded(ttt);
    Ap                  = vAPathExpanded(ttt+1);
    

    % B. Key values
    [~,C,~,~]           = fnFindNTransitions(K,Kp,A,Parameters);
    [~,Cp,Rp,~]       = fnFindNTransitions(Kp,Kpp,Ap,Parameters);

    % C. Euler error 
    RHS                 = C^(-pSigma) * (1 + Psi_2(K,Kp));
    LHS                 = pBeta * Cp^(-pSigma) * (1 + Rp - Psi_1(Kp,Kpp));
    vResiduals(ttt)     = RHS - LHS;
end 

end 