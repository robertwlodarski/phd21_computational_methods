% Calibrate using the pattern search algorithm

% Solve for
% 1. pChi
% 2. pEta
% 3. pb
% 4. pSigmaz

% Function
tic;
vInitial                    = [0.956,4.608,0.062,0.276];
fnMinimise                  = @(x) fnLossFunctionQuadratic(vTargetMoments,mWeights,pTau,x(2),x(1),pBeta,x(4),pAlpha,pA,pa,pr,x(3),pMaxIter,pStepSize);
fnOptions                   = optimset('Display','iter','MaxTime',3*3600);
[vParamsSolvedPS, LossPS]   = patternsearch(fnMinimise,vInitial,[],[],[],[],vLowerBound,vUpperBound,[],fnOptions);
ElapsedTimePS               = toc / 60;         

% Display time message
fprintf('Elapsed time for standard calibration: %.2f minutes\n', ElapsedTimePS)

% Parameters
pChiPS                      = vParamsSolvedPS(1);
pEtaPS                      = vParamsSolvedPS(2);
pbPS                        = vParamsSolvedPS(3);
pSigmazPS                   = vParamsSolvedPS(4);

% Save results 
save('_results/ParametersPS.mat', 'pChiPS', 'pEtaPS', 'pbPS', 'pSigmazPS', 'LossPS');