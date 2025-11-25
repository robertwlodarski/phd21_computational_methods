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

%% 2. Run the loop
for ttt = 1:1:pT
    
    % A. K, Kp, Kpp, A, and Ap
    % t=1 case
    if ttt == 1
        K               = K0; 
        Kp              = vKGuess(1);
        Kpp             = vKGuess(2);
        A               = vAPath(1);
        Ap              = vAPath(2);
    elseif ttt == pT
        K               = vKGuess(ttt-1);
        Kp              = vKGuess(ttt);
        Kpp             = KT;
        A               = vAPath(ttt);
        Ap              = vAPath(end);
    else 
        K               = vKGuess(ttt-1);
        Kp              = vKGuess(ttt);
        Kpp             = vKGuess(ttt+1);
        A               = vAPath(ttt);
        Ap              = vAPath(ttt+1);
    end 
    
    % B. Key values
    [~,C,~,~]           = fnFindNTransitions(K,Kp,A,Parameters);
    [~,Cp,Rp,~]         = fnFindNTransitions(Kp,Kpp,Ap,Parameters);

    % C. Euler error 
    RHS                 = C^(-pSigma) * (1 + Psi_2(K,Kp));
    LHS                 = pBeta * Cp^(-pSigma) * (1 + Rp - Psi_1(Kp,Kpp));
    vResiduals(ttt)     = LHS - RHS;
end 

end 

