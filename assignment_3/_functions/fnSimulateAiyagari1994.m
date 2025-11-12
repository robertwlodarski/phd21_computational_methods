function [RESULT]       = fnSimulateAiyagari1994(Parameters,Grids)

%% 1. Unpacking
% Parameters
pAlpha                  = Parameters.pAlpha;
pDelta                  = Parameters.pDelta;
pA                      = Parameters.pA;
% Grids
vGridA1                 = Grids.vGridA1;
vGridZ                  = Grids.vGridZ;
mTransitionZ            = Grids.mTransitionZ;


%% 2. Set-up
% Iterations business: Convergence, acceleration, and GE
iWeightOld              = 0.9;
iErrorGE                = 10;
iTolGE                  = 1e-8;
iTolVFI                 = 1e-8;
iAccelerationInterv     = 20;
iAccelerationStart      = 30;
iIterNumGE              = 1;
% Matrices: Value function, policy functions
iValue                  = repmat(0.1*vGridA1,1,size(vGridZ,1));
iValueNew               = zeros(size(iValue));
iPolicyCons             = iValueNew;
iPolicyWealthNext       = iPolicyCons;
iExpectedValue          = iValue * mTransitionZ';
% Distributions
iCurrentDistribution    = ones(size(iValue)) / sum(sum(ones(size(iValue))));

%% 3. Stationary dist. of labour allocation
iErrorZDist             = 10;
iTolZDist               = 1e-8;
vZDist                  = ones(size(vGridZ,1),1) / sum( ones(size(vGridZ,1),1));
while iErrorZDist>iTolZDist
    vZDistNext          = mTransitionZ' * vZDist;
    iErrorZDist         = abs(vZDistNext-vZDist);
    vZDist              = vZDistNext; 
end
vLabSupply              = vGridZ'*vZDist;

%% 4. GE & VFI loops

% Initial K guess
iK                      = 6;
% START GE LOOP
while iErrorGE>iTolGE
    % Derive K-related items
    iInterest           = pAlpha * pA * (iK / vLabSupply)^(pAlpha - 1) - pDelta;
    iWage               = (1 - pAlpha) * pA * (iK / vLabSupply)^(pAlpha);
    % Prepare for starting VFI loop
    iErrorVFI           = 10;
    iNumIterVFI         = 1;
    % START VFI LOOP
    while iErrorVFI>iTolVFI
        % Start iterating over z and a
        for zzz = 1:1:size(vGridZ,1)
            % Set z-related items
            z           = vGridZ(zzz);
            % Expected future values
            iExpValueZ  = iExpectedValue(:,zzz);
            iMinWealth  = vGridA1(1);
            for aaa = 1:1:size(vGridA1,1) %FINISHED HERE
            end
        end
        % Update VFI loop items
        iNumIterVFI     = iNumIterVFI + 1;
        
        % END VFI LOOP
    end
    % END GE LOOP
end

end