function [Res]           = fnEulerResidual(K,Kp,vK,vKNextGuess,EulerError,BCError,ConsImplied,Interest,Wage)

% A. Find contemporaneous N, W, C and R
[N,~,~,~]               = fnFindN(K,Kp,BCError,ConsImplied,Interest,Wage);

% B. Fast interpolation of K''
Kpp_bot                        = sum(vK<=Kp);
Kpp_bot(Kpp_bot<1)             = 1;
Kpp_bot(Kpp_bot>=size(vK,1))   = size(vK,1)-1;
Kpp_up                         = Kpp_bot+1;
Bot_weight                     = (vK(Kpp_up)-Kp)/(vK(Kpp_up)-vK(Kpp_bot));
Kpp                            = Bot_weight * vKNextGuess(Kpp_bot) + (1 - Bot_weight) * vKNextGuess(Kpp_up);

% C. Get tomorrow's elements
[Np,~,~,~]              = fnFindN(Kp,Kpp,BCError,ConsImplied,Interest,Wage);

% D. Get the residual
Res                     = EulerError(N,K,Kp,Np,Kpp);

end 