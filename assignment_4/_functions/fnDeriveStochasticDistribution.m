function [Distribution] = fnDeriveStochasticDistribution(vAssets,vLabour,pN,Grids)
    % Load grids 
    vGridA1             = Grids.vGridA1;
    vGridZ              = Grids.vGridZ;

    % Iterate
    Distribution        = zeros(size(vGridA1,1),size(vGridZ,1));
    
    % Start iterating over z and a
    for nnn = 1:1:size(vAssets,1)
        % Set up
        iWealth                         = vAssets(nnn);
        izIndex                         = vLabour(nnn);
        iLB                             = sum(vGridA1<iWealth);
        iLB(iLB<=0)                     = 1;
        iLB(iLB>=size(vGridA1,1))       = size(vGridA1,1)-1;
        iUB                             = iLB + 1;
        iWeightLB                       = (vGridA1(iUB)-iWealth) / (vGridA1(iUB)-vGridA1(iLB));
        iWeightLB(iWeightLB<0)          = 0;
        iWeightLB(iWeightLB>1)          = 1;
        iWeightUB                       = 1 - iWeightLB;
        iMass                           = 1;

        % Upper
        Distribution(iUB,izIndex)  = Distribution(iUB,izIndex) + iMass * iWeightUB;
        % Lower
        Distribution(iLB,izIndex)  = Distribution(iLB,izIndex) + iMass * iWeightLB;
    end
Distribution            = Distribution / pN;
end