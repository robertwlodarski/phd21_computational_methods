function [w, T] = fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize)
    % Initialise
    iWage           = 1.0;
    iT              = 1.0;
    iIteration      = 0;
    iConvergence    = false;

    % Start the process
    while iIteration < pMaxIter
        % Households supply
        [iEffectiveLab, iEmp, ~]    = fnAggregateLabourSupply(iWage,iT,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize);
        iWage                       = (iEffectiveLab<1e-8) * 2 * iWage + (iEffectiveLab>1e-8)*iWage;
        
        % Firm demand & transfers
        iWageNext           = fnLabourDemand(pa,iEffectiveLab,pA,pAlpha);
        iTNext              = fnGovernmentLumpTransfer(iWage,iEffectiveLab, iEmp,pa,pr,pb,pTau);

        % Update error
        iDeviation          = [iWage - iWageNext; iT - iTNext];
        iError              = iDeviation' * iDeviation;

        % Convergence criterion
        if iError < 1e-5
            iConvergence    = true;
            break
        else 
            % If no convergence, update the initial w & T
            iWage           = pStepSize * iWageNext + (1 - pStepSize) * iWage;
            iT              = pStepSize * iTNext + (1 - pStepSize) * iT;        
        end
        iIteration  = iIteration + 1;
        
    % End the loop
    end
    
    % Save the results
    w               = iWage;
    T               = iT;

    % Print message(s)
    if ~iConvergence
    fprintf('STATUS: No satisfying convergence. \n')
    else 
    fprintf('STATUS: Solution found. \n')
    end 
    fprintf('Quadratic error    = %.5f\n', iError)
    fprintf('# of iterations    = %.0f\n', iIteration)
    fprintf('Wage               = %.2f\n', w)
    fprintf('Transfer           = %.2f\n',T)
end