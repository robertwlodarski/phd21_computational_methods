function [vGrid, mTransition] = fnTauchenLogNormalZ(Parameters)
    % Log-normal structure 
    % X = log(Z)
    % Compute grid for X, and then use Z=exp(X)
    % P remains the same

    % Unpack parameters
    N               = Parameters.pTauchenN;
    M               = Parameters.pTauchenM;
    SigmaZ          = Parameters.pSigmaZ;
    RhoZ            = Parameters.pRhoZ;

    % Grid preparation
    vXGrid = zeros(M,1);
    vXGrid(1)       = -N * SigmaZ / sqrt(1 - RhoZ^2);
    vXGrid(end)     = N * SigmaZ / sqrt(1 - RhoZ^2);
    d               = (vXGrid(end)-vXGrid(1)) / (M-1);
    for iii         = 2:size(vXGrid,1)
        vXGrid(iii) = vXGrid(iii-1)+d;
    end
    vGrid           = exp(vXGrid);    

    % Transition probabilities 
    mTransition     = zeros(M,M);
    for jjj = 1:1:M
        for iii = 1:1:M
            IntPlus                 = (vXGrid(jjj) + d/2 - RhoZ *vXGrid(iii)) / SigmaZ;
            IntMinus                = (vXGrid(jjj) - d/2 - RhoZ *vXGrid(iii)) / SigmaZ;
            P                       = normcdf(IntPlus,0,1) -  normcdf(IntMinus,0,1);
            mTransition(iii,jjj)    = P;
        end
    end
    
    % Ensure transition rows sum to 1 (they don't for N=1)
    mTransition             = mTransition ./ sum(mTransition,2);
end 