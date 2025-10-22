function T = fnGovernmentLumpTransfer(w,eff_labour_supply, employment,pa,pr,pb,pTau)
    T = (eff_labour_supply * w + pa * pr) * pTau - (1 - employment) * pb;
end