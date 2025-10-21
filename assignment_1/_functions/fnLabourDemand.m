function w = fnLabourDemand(K,L,pA,pAlpha)
    w           = pA * (1 - pAlpha) * (K / L)^pAlpha;
end