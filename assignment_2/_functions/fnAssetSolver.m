function [a_prime, c]       = fnAssetSolver(a,h,z,t,working,mV,mS,Parameters, Grid)
    
    % 1. Unpack parameters
    vGridZ                  = Grids.vGridZ;
    vGridH                  = Grids.vGridH;
    vGrida                  = Grids.vGrida;
    vGridAge                = Grids.vGridAge;    
    
    % 2. Set interpolants 
    pInterMethod            = {'linear', 'nearest', 'nearest', 'nearest', 'linear'};
    vInterGrids             = {vGrida,vGridH,vGridZ,vGridAge,vGrida};
    F_mV                    = griddedInterpolant(vInterGrids,mV);
 
    % 3. Design function to be maximised
    if working == 1
    else
    end 

end

