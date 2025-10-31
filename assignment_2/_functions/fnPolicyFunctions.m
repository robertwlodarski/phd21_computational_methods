function [mWorkingPolicyMatrix, mNonWorkingPolicyMatrix] = fnPolicyFunctions(Parameters,Grids)
    % Unpack the key elements
    pEta                    = Parameters.pEta;
    pBeta                   = Parameters.pBeta;
    pr                      = Parameters.pr;
    pb                      = Parameters.pb;
    pVarphi                 = Parameters.pVarphi;
    pPhi                    = Parameters.pPhi;
    vGridZ                  = Grids.vGridZ;
    vGridH                  = Grids.vGridH;
    vGrida                  = Grids.vGrida;
    vGridAge                = Grids.vGridAge;

    % Initialise policy function matrices
    % (a,h,z,T) -> a'
    mWorkingPolicyMatrix    = zeros(size(vGrida,1),size(vGridH,1),size(vGridZ,1),size(vGridAge,1)); 
    mNonWorkingPolicyMatrix = mWorkingPolicyMatrix;

    % Initialise value function matrices
    % (a,h,z,T)
    mV                      = zeros(size(vGrida,1),size(vGridH,1),size(vGridZ,1),size(vGridAge,1));
    mS                      = mV;
    mW                      = mS;
    mN                      = mW;

    % Solve for T=50
    for aaa = 1:1:size(vGrida,1)
        for hhh = 1:1:size(vGridH)
            for zzz = 1:1:size(vGridZ)
                % Wages, assets and human capital
                h           = vGridH(hhh);
                a           = vGrida(aaa);
                z           = vGridZ(zzz);
                w           = fnWage(h,z,Parameters);
                % Working household
                a_next_w    = (pBeta/(pBeta+1))*(w + (1+pr)*a);
                c_w         = (1 / pBeta) * a_next_w;
                W           = log(c_w) - pEta + pBeta * log(a_next_w);
                % Non-working household
                a_next_n    = (pBeta/(pBeta+1))*(pb + (1+pr)*a);
                c_n         = (1 / pBeta) * a_next_n;
                N           = log(c_n) + pBeta * log(a_next_n);
                % Value functions
                V           = fnGumbelTrick(W,N,Parameters);
                arg1        = pVarphi * W + (1 - pVarphi) * N - pPhi;
                S           = fnGumbelTrick(arg1,N,Parameters);

            end
        end 
    end



end