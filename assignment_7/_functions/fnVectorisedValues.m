function    [vV, vC, vN]        = fnVectorisedValues(vAp,vA,iLabourProd,iWage,iInterest,Parameters,Grids,iExpValueZ)

%% 1. Initialisation
vGridA1             = Grids.vGridA1;
pBeta               = Parameters.pBeta;
pMu                 = Parameters.pMu;
pEta                = Parameters.pEta;
pChi                = Parameters.pChi;

%% 2. Approximation
% A. Discretise 
iIndices                        = discretize(vAp,[-Inf; vGridA1; Inf]);

% B. Top & bottom indices
iIndLow                             = iIndices-1;
iIndLow(iIndLow <1)                 = 1;
iIndLow(iIndLow >= size(vGridA1,1)) = size(vGridA1,1)-1;
iIndHigh                            = iIndLow+1;
iWeightLow                          = (vGridA1(iIndHigh) - vAp) ./ (vGridA1(iIndHigh)-vGridA1(iIndLow));
iWeightLow(iWeightLow < 0)          = 0;
iWeightLow(iWeightLow > 1)          = 1;

% C. Expected values 
vEVLow                              = iExpValueZ(iIndLow);
vEVHigh                             = iExpValueZ(iIndHigh);
vEV                                 = iWeightLow .* vEVLow + (1 - iWeightLow) .* vEVHigh;

%% 3. Produce final results: Compute V
A                                   = (1 + iInterest) * vA - vAp - pMu / 2 * (vAp.^2 ./ vA - 2*vAp + vA);
B                                   = (iWage * iLabourProd)^2 / pEta;
vC                                  = (sqrt(A.^2+ 4 * B)+A) / 2;
vN                                  = iWage * iLabourProd ./ (vC * pEta);
vV                                  = log(vC) - pEta * vN .^(1 + 1 / pChi) / (1 + 1 / pChi) + pBeta * vEV; 

end 