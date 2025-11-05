function [vGrid, mTransition] = fnTauchenLogNormal(Parameters)
    % Log-normal structure 
    % X = log(Z)
    % Compute grid for X, and then use Z=exp(X)
    % P remains the same

    % Unpack parameters
    N               = Parameters.pTauchenN;
    M               = Parameters.pTauchenM;
    Sigmaz          = Parameters.pSigmaz;

    % Grid preparation
    vXGrid = zeros(M,1);
    vXGrid(1)       = -N * Sigmaz;
    vXGrid(end)     = N * Sigmaz;
    d               = (vXGrid(end)-vXGrid(1)) / (M-1);
    for iii         = 2:size(vXGrid,1)
        vXGrid(iii) = vXGrid(iii-1)+d;
    end
    vGrid           = exp(vXGrid);    

    % Transition probabilities 
    mTransition         = zeros(M,M);
    for jjj = 1:1:M
        IntPlus             = (vXGrid(jjj) + d/2) / Sigmaz;
        IntMinus            = (vXGrid(jjj) - d/2) / Sigmaz;
        P                   = normcdf(IntPlus,0,1) -  normcdf(IntMinus,0,1);
        mTransition(:,jjj)  = repmat(P,[M,1]);
    end
    
end 