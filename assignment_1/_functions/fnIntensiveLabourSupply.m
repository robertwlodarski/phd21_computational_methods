function [c,n]      = fnIntensiveLabourSupply(w,T,pz,pTau,pEta,pChi,pBeta,pa,pr,pMaxIter,pStepSize)
    c_init          = 1.0;          % Initial guess of c
    c_temp          = c_init;       % Initialise c_temp
    conv            = false;
    iter            = 0;

    % Start iterations 
    while iter < pMaxIter
        RHS         = fnIntensiveLabourSupplyRHS(w,c_temp,T,pz,pTau,pEta,pChi,pa,pr);
        % Reduce c_temp if the RHS value is too large
        if RHS < 0
            c_temp  = 0.5 * c_temp;
            continue
        end 
        c_new       = RHS / (1+pBeta);
        err         = abs(c_new-c_temp);
        % Stop if the error is small
        if err < 1e-5
            conv    = true;
            break
        else 
            % Update the iteration
            c_temp  = pStepSize * c_new + (1 - pStepSize) * c_temp;  
        end 
        iter        = iter + 1;
    end
    % Print out error
    if ~conv
        error("Euler equation did not converge.")
    end 
    % Save the final values
    c               = c_temp;
    n               = ((pz * w * (1 - pTau))/( c * pEta))^pChi;
end 