function [m1, m2, m3, m4] = fnComputeMoments(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize)
    % Find the equilibrium
    [w, T]                                                  = fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize);
    [eff_labour_supply, employment, eff_labour_supply_std]  = fnAggregateLabourSupply(w,T,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize);
    unemployment                                            = 1 - employment;

    % Derive the moments
    % m1            = total working hours
    % m2            = unemployment rate
    % m3            = benefit as per cent of wage
    % m4            = STD/ Mean of labour income
    m1              = eff_labour_supply / employment;
    m2              = unemployment;
    m3              = pb / (eff_labour_supply * w);
    m4              = eff_labour_supply_std / eff_labour_supply;
end