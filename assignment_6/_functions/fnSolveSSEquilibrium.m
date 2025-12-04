function [Results]      = fnSolveSSEquilibrium(Parameters)
    % Unpack parameters 
    pAlpha              = Parameters.pAlpha;
    pDelta              = Parameters.pDelta;
    pBeta               = Parameters.pBeta;
    pEta                = Parameters.pEta;
    pChi                = Parameters.pChi;
    
    % Constants
    R                   = 1 / pBeta -1;
    K_N                 = (pAlpha / (R + pDelta))^(1 / (1-pAlpha));
    W                   = (1-pAlpha) * K_N^(pAlpha);
    C_N                 = W + R * K_N;
    N                   = (W / (C_N * pEta))^(pChi / (1+pChi));
    C                   = C_N * N;
    K                   = K_N *N;

    % Save results
    Results             = struct('R', R, 'K', K, 'W', W, 'C', C, 'N', N);
end 