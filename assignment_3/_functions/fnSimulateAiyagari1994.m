function [RESULT]       = fnSimulateAiyagari1994(Parameters,Grids)

%% 1. Unpacking
% Parameters

% Grids
vGridA1                 = Grids.vGridA1;
vGridZ                  = Grids.vGridZ;


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
% Distributions
iCurrentDistribution    = ones(size(iValue)) / sum(sum(ones(size(iValue))));


% START GE LOOP
while iErrorGE>iTolGE
    % Derive K-related items
    % Prepare for starting VFI loop
    iErrorVFI           = 10;
    iNumIterVFI         = 1;
    % START VFI LOOP
    while iErrorVFI>iTolVFI
        
        % Update VFI loop items
        iNumIterVFI     = iNumIterVFI + 1;
        
        % END VFI LOOP
    end
    % END GE LOOP
end

end