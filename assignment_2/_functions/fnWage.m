function w          = fnWage(h,z,Parameters)
    % Unpack
    pGammah         = Parameters.pGammah;
    pGammaz         = Parameters.pGammaz;
    pGamma0         = Parameters.pGamma0;

    % Compute the wage
    w = pGamma0 + pGammah * h + pGammaz * z;
end 