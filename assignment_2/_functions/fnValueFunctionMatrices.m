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
                a_next_w    = max(a_next_w, 1e-9);
                c_w         = (1 / pBeta) * a_next_w;
                c_w         = max(c_w, 1e-9);
                W           = log(c_w) - pEta + pBeta * log(a_next_w);
                % Non-working household
                a_next_n    = (pBeta/(pBeta+1))*(pb + (1+pr)*a);
                a_next_n    = max(a_next_n,1e-9);
                c_n         = (1 / pBeta) * a_next_n;
                c_n         = max(c_n,1e-9);
                N           = log(c_n) + pBeta * log(a_next_n); 
                % Value functions
                V           = fnGumbelTrick(W,N,Parameters);
                arg1        = pVarphi * W + (1 - pVarphi) * N - pPhi;
                S           = fnGumbelTrick(arg1,N,Parameters);
                % Save the "incorrect choice" results
                mV(aaa,hhh,zzz,pT+1,:) = -Inf;
                mS(aaa,hhh,zzz,pT+1,:) = -Inf;
                mW(aaa,hhh,zzz,pT+1,:) = -Inf;
                mN(aaa,hhh,zzz,pT+1,:) = -Inf;
                % Find the index of a' that satisies the correct result
                [~ , index_ap_w]       = min(abs(vGrida - a_next_w));
                [~ , index_ap_n]       = min(abs(vGrida - a_next_n));
                % Assign better values
                mV(aaa,hhh,zzz,pT+1,index_ap_w) = V;
                mS(aaa,hhh,zzz,pT+1,index_ap_n) = S;
                mW(aaa,hhh,zzz,pT+1,index_ap_w) = W;
                mN(aaa,hhh,zzz,pT+1,index_ap_n) = N;
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
                        h           = vGridH(hhh);
                        a           = vGrida(aaa);
                        z           = vGridZ(zzz);
                        w           = fnWage(h,z,Parameters);                        
                        % Working household
                        c_w         = w + (1+pr)*a - a_next;
                        c_w         = max(c_w, 1e-9);
                        W           = log(c_w) - pEta + pBeta * (mTransition(zzz,:) * max(squeeze(mV(aap,h_fut,:,ttt+1,:)),[],2));
                        % Non-working household
                        c_n         = pb + (1+pr)*a - a_next;
                        c_n         = max(c_n,1e-9);
                        N           = log(c_n) + pBeta * (mTransition(zzz,:) * max(squeeze(mS(aap,hhh,:,ttt+1,:)),[],2)); 
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