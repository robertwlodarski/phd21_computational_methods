function Results    = fnSteadyState(Parameters)

%% 1. Unpack parameters to be used
pBeta               = Parameters.pBeta;
pDelta              = Parameters.pDelta;
pAlpha              = Parameters.pAlpha;
pEta                = Parameters.pEta;
pGamma              = Parameters.pGamma;

%% 2. Compute important values
vK_to_Y             = pAlpha / (1 / pBeta - (1 - pDelta));
vN                  = pGamma / (pEta * (1 - pDelta * vK_to_Y));
vK                  = (vK_to_Y * vN^(pGamma))^(1 / (1-pAlpha));
vY                  = vK^pAlpha * vN^pGamma;
vW                  = pGamma * vY / vN;
vC                  = 1 / pEta * vW;

%% 3. Save the results
Results.N           = vN;
Results.K           = vK; 
Results.Y           = vY; 
Results.C           = vC; 
Results.W           = vW;

end 