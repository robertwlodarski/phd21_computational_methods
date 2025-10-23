function ExcessSupply      =    fnExcessLabourSupply(w,T,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize,pA,pAlpha)
    % Labour supply & demand
    [Supply,~,~]            = fnAggregateLabourSupply(w,T,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize);
    Demand                  = fnLabourDemandInverse(pa,w,pA,pAlpha);
    % Excess value
    ExcessSupply            = Supply - Demand;
end
