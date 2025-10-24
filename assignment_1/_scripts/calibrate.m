% Model calibration

% Solve for
% 1. pChi
% 2. pEta
% 3. pb
% 4. pSigmaz

% Function
tic;
vInitial                = log([0.7744,6.4986,0.0633,0.3136]);
fnMinimise              = @(x) fnLossFunctionQuadratic(vTargetMoments,mWeights,pTau,exp(x(2)),exp(x(1)),pBeta,exp(x(4)),pAlpha,pA,pa,pr,exp(x(3)),pMaxIter,pStepSize);
fnOptions               = optimset('Display','iter');
[vParamsSolved, Loss]   = fminsearch(fnMinimise,vInitial,fnOptions);
ElapsedTime             = toc / 60;   

% Display time message
fprintf('Elapsed time for standard calibration: %.2f minutes\n', ElapsedTime);

% Parameters
pChiBas                 = exp(vParamsSolved(1));
pEtaBas                 = exp(vParamsSolved(2));
pbBas                   = exp(vParamsSolved(3));
pSigmazBas              = exp(vParamsSolved(4));
LossBas                 = Loss;

% Parameters: (27 minutes)
% Chi        = 0.956 
% b          = 0.062 
% Eta        = 4.608 
% Sigma_z    = 0.276 
% Q. error   = 0.11310 

save('_results/ParametersBas.mat', 'pChiBas', 'pEtaBas', 'pbBas', 'pSigmazBas', 'LossBas');