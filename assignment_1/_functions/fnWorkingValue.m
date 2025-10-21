function W = fnWorkingValue(c_eq,n_eq,pBeta,pChi,pEta)
    % Calculate assets & Frisch elasticity term
    a_prime             = pBeta * c_eq;
    frisch              = 1 + 1/pChi;
    % Obtain the value function
    W                   = log(c_eq) - pEta * (1 / frisch) * n_eq^frisch +  pBeta * log(a_prime);
end 