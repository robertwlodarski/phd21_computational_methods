function [N,C,R,W]  = fnFindN(K,Kp,BCError,ConsImplied,Interest,Wage)
    % 1. Solve for N
    BCErrorKSq      = @(N) BCError(N,K,Kp)^2;
    N               = fminbnd(BCErrorKSq,1e-5,3-1e-5);

    % 2. Solve for other elements
    C               = ConsImplied(N,K);
    R               = Interest(N,K);
    W               = Wage(N,K);
 
end 