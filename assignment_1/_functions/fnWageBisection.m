function w_star         = fnWageBisection(T,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize,pA,pAlpha,pMaxWage)
    % Excess supply
    fnExcessSupply      = @(w) fnExcessLabourSupply(w,T,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize,pA,pAlpha);
    fnExcessSupplyVec   = @(w) arrayfun(fnExcessSupply,w);
    % Grid
    vGrid               = linspace(0.01,pMaxWage,1000);
    vExcessSupplyVals   = fnExcessSupplyVec(vGrid);
    % Threshold values
    vLower              = sum((vExcessSupplyVals<0));
    vUpper              = vLower+1;
    % Calculation
    vWeight             = vExcessSupplyVals(vUpper) / (vExcessSupplyVals(vUpper)-vExcessSupplyVals(vLower));
    w_star              = vGrid(vLower) * vWeight + vGrid(vUpper) * (1 - vWeight);
end 