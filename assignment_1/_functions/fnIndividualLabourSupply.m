function h = fnIndividualLabourSupply(w,T,pTau,pEta,pChi,pBeta,pa,pr,pz,pb,pMaxIter,pStepSize)
    % "If working" consumption and hours 
    [c_iw, n_if]            = fnIntensiveLabourSupply(w,T,pz,pTau,pEta,pChi,pBeta,pa,pr,pMaxIter,pStepSize);
    % Value functions
    WW                      = fnWorkingValue(c_iw,n_if,pBeta,pChi,pEta);
    NN                      = fnNonWorkingValue(T,pb,pa, pr, pTau, pBeta);
    % Value of working 
    if WW >= NN
        h = n_if;
    else 
        h = 0.0;
    end
end