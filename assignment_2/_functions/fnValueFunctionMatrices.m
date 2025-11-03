function [mV,mS,mW,mN]      = fnValueFunctionMatrices(Parameters,Grids)
    % Unpack the key elements
    pEta                    = Parameters.pEta;
    pBeta                   = Parameters.pBeta;
    pr                      = Parameters.pr;
    pb                      = Parameters.pb;
    pVarphi                 = Parameters.pVarphi;
    pPhi                    = Parameters.pPhi;
    pT                      = Parameters.pT;
    vGridZ                  = Grids.vGridZ;
    vGridH                  = Grids.vGridH;
    vGrida                  = Grids.vGrida;
    vGridAge                = Grids.vGridAge;
    mTransition             = Grids.mTransitionZ;

    % Initialise value function matrices
    % (a,h,z,T,a')
    mV                      = zeros(size(vGrida,1),size(vGridH,1),size(vGridZ,1),size(vGridAge,1),size(vGrida,1));
    mS                      = mV;
    mW                      = mS;
    mN                      = mW;

    % Solve for T=50
    for aaa = 1:1:size(vGrida,1)
        for hhh = 1:1:size(vGridH,1)
            for zzz = 1:1:size(vGridZ,1)
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
                % Save the results
                mV(aaa,hhh,zzz,pT+1,:) = V;
                mS(aaa,hhh,zzz,pT+1,:) = S;
                mW(aaa,hhh,zzz,pT+1,:) = W;
                mN(aaa,hhh,zzz,pT+1,:) = N;
            end
        end 
    end

    % Solve for the remaining times (backward)
    for ttt=pT:(-1):1                               % Time
        for aaa = 1:1:size(vGrida,1)                % Existing wealth
            for hhh = 1:1:size(vGridH,1)            % Human capital
                for zzz = 1:1:size(vGridZ,1)        % Productivity
                    for aap = 1:1:size(vGrida,1)    % Future wealth
                        % Wages, assets and human capital
                        a_next      = vGrida(aap);
                        h_fut       = min(size(vGridH,1),hhh+1);
                        % Working household
                        c           = (1 / pBeta) * a_next;
                        W           = log(c) - pEta + pBeta * (mTransition(zzz,:) * squeeze(mV(aaa,h_fut,:,ttt+1,aap)));
                        % Non-working household
                        N           = log(c) + pBeta * (mTransition(zzz,:) * squeeze(mS(aaa,h_fut,:,ttt+1,aap)));
                        % Value functions
                        V           = fnGumbelTrick(W,N,Parameters);
                        arg1        = pVarphi * W + (1 - pVarphi) * N - pPhi;
                        S           = fnGumbelTrick(arg1,N,Parameters);
                        % Save the results
                        mV(aaa,hhh,zzz,ttt,aap) = V;
                        mS(aaa,hhh,zzz,ttt,aap) = S;
                        mW(aaa,hhh,zzz,ttt,aap) = W;
                        mN(aaa,hhh,zzz,ttt,aap) = N;
                    end
                end 
            end 
        end 
    end

end