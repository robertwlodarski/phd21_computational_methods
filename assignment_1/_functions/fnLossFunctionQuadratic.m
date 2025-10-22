function Loss           = fnLossFunctionQuadratic(vTargetMoments,mWeights,pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize)
    try
        % Compute model moments
        [m1, m2, m3, m4] = fnComputeMoments(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize);
        vModelMoments   = [m1; m2; m3; m4];
        % Deviation & loss
        vDeviation      = (vModelMoments - vTargetMoments) ./ vTargetMoments;
        Loss            = vDeviation' * mWeights * vDeviation;
        % Printing results
        fprintf('Parameters: \n')
        fprintf('Chi        = %.3f \n', pChi)
        fprintf('b          = %.3f \n', pb)
        fprintf('Eta        = %.3f \n', pEta)
        fprintf('Sigma_z    = %.3f \n', pSigmaz)
        fprintf('Q. error   = %.5f \n', Loss)
    catch 
        % Print that there was an issue
        fprintf('Error occurred during computation.\n');
        Loss            = 1e4;
    end
end