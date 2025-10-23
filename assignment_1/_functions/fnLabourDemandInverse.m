function L = fnLabourDemandInverse(K,w,pA,pAlpha)
    L           = ((1-pAlpha) * pA / w)^(1 / pAlpha) * K;
end