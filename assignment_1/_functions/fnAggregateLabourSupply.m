function [eff_labour_supply, employment, eff_labour_supply_std] = fnAggregateLabourSupply(w,T,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize)
    % Set upper and lower values for z (0 and inverse CDF)
    z_lower                 = 0;
    z_upper                 = logninv(1-1e-10);
    % Labour functions
    fn_Extensive            = @(z) fnExtensiveLabourSupply(w,T,pTau,pEta,pChi,pBeta,pa,pr,z,pb,pMaxIter,pStepSize);
    fn_Individual           = @(z) fnIndividualLabourSupply(w,T,pTau,pEta,pChi,pBeta,pa,pr,z,pb,pMaxIter,pStepSize);
    % Allow for working with vectors
    fn_Extensive_vec        = @(z) arrayfun(fn_Extensive,z);
    fn_Individual_vec       = @(z) arrayfun(fn_Individual,z);
    % Effective labour supply
    fn_eff_integrand        = @(z) lognpdf(z,0,pSigmaz) .* z .* fn_Individual_vec(z);
    eff_labour_supply       = integral(fn_eff_integrand,z_lower,z_upper);
    % Employment 
    fn_emp_integrand        = @(z)  lognpdf(z,0,pSigmaz) .* fn_Extensive_vec(z);
    employment              = integral(fn_emp_integrand, z_lower, z_upper);
    % Effective labour supply standard deviation
    fn_eff_std_integrand    = @(z) lognpdf(z,0,pSigmaz) .* (z .* fn_Individual_vec(z)).^2;
    eff_labour_supply_2m    = integral(fn_eff_std_integrand, z_lower, z_upper);
    eff_labour_supply_std   = sqrt(eff_labour_supply_2m - eff_labour_supply^2);
end 