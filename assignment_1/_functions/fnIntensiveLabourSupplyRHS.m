function [RHS] = fnIntensiveLabourSupplyRHS(w,c,T,pz,pTau,pEta,pChi,pa,pr)
    inner_term      = pz * w * (1 - pTau);
    RHS             = inner_term * (inner_term / (pEta * c))^(pChi) + pa * (1 + pr * (1 - pTau)) + T;
end 