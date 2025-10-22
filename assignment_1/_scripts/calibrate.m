% Model calibration

% Solve for
% 1. pChi
% 2. pEta
% 3. pb
% 4. pSigmaz

% Function
tic;
vInitial                = [0,log(1.1507),log(0.1230),log(0.2996)];
fnMinimise              = @(x) fnLossFunctionQuadratic(vTargetMoments,mWeights,pTau,exp(x(2)),exp(x(1)),pBeta,exp(x(4)),pAlpha,pA,pa,pr,exp(x(3)),pMaxIter,pStepSize);
fnOptions               = optimset('Display','iter');
[vParamsSolved, Loss]   = fminsearch(fnMinimise,vInitial,fnOptions);
ElapsedTime             = toc / 60;         

% Parameters
pChi                    = exp(vParamsSolved(1));
pEta                    = exp(vParamsSolved(2));
pb                      = exp(vParamsSolved(3));
pSigmaz                 = exp(vParamsSolved(4));

% Parameters: (27 minutes)
% Chi        = 0.956 
% b          = 0.062 
% Eta        = 4.608 
% Sigma_z    = 0.276 
% Q. error   = 0.11310 