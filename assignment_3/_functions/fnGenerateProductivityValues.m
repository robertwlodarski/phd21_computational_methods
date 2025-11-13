function [mProductivityValues]      = fnGenerateProductivityValues(mProductivityShocks, Grids)
    % Unpacking 
    mTransitionZ                        = Grids.mTransitionZ;         
    mTransitionZCDF                     = cumsum(mTransitionZ,2);

    % Find the marginal distribution
    v0                                  = ones(size(mTransitionZ,1),1)./size(mTransitionZ,1);
    vMarginal                           = (mTransitionZ')^100 *  v0;
    vMarginalCDF                        = cumsum(vMarginal);

    % First period 
    mProductivityValues                 = zeros(size(mProductivityShocks));
    mProductivityValues(1,:)            = sum(mProductivityShocks(1,:) > vMarginalCDF)+1;

    % Remaining periods
    for iii = 2:1:size(mProductivityValues,1)
        for jjj = 1:1:size(mProductivityValues,2)
            z_prev                      = mProductivityValues(iii-1,jjj);
            shock_v                     = mProductivityShocks(iii,jjj);
            index                       = sum(shock_v > mTransitionZCDF(z_prev,:))+1;
            mProductivityValues(iii,jjj)=index;
        end
    end
end