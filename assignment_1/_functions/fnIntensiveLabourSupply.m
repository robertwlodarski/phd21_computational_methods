function [c,n]      = fnIntensiveLabourSupply(w,T,pz,pTau,pEta,pChi,pBeta,pa,pr,pMaxIter,pStepSize)
    iConsumption    = 1.0;       % Initialise consumpton
    iConvergence    = false;
    iIteration      = 0;

    % Start iterations 
    while iIteration < pMaxIter
        iRHS        = fnIntensiveLabourSupplyRHS(w,iConsumption,T,pz,pTau,pEta,pChi,pa,pr);
        % Reduce c_temp if the RHS value is too large
        if iRHS < 0
            iConsumption  = 0.5 * iConsumption;
            continue
        end 
        iConsumptionNext    = iRHS / (1+pBeta);
        iError                 = abs(iConsumptionNext-iConsumption);
        % Stop if the error is small
        if iError < 1e-5
            iConvergence    = true;
            break
        else 
            % Update the iteration
            iConsumption  = pStepSize * iConsumptionNext + (1 - pStepSize) * iConsumption;  
        end 
        iIteration        = iIteration + 1;
    end
    % Print out error
    if ~iConvergence
        error("Euler equation did not converge.")
    end 
    % Save the final values
    c               = iConsumption;
    n               = ((pz * w * (1 - pTau))/( c * pEta))^pChi;
end 