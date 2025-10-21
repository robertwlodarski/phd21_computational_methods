function NW = fnNonWorkingValue(T,pb,pa, pr, pTau, pBeta)
    % Get the consumption from the non-working Euler
    numerator           = pb + pa * (1 + pr * (1 - pTau)) + T;
    c_eq                = numerator / (1+pBeta);
    % Tomorrow's assets
    a_prime             = pBeta * c_eq;
    % Value function
    NW                  = log(c_eq) +  pBeta * log(a_prime);    
end