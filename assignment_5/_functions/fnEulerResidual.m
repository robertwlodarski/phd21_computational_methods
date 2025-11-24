function [Res]          = fnEulerResidual(K,Kp,vK,vKNextGuess,EulerError,BCError,ConsImplied,Interest,Wage)

% A. Find contemporaneous N, W, C and R
[N,C,~,~]               = fnFindN(K,Kp,BCError,ConsImplied,Interest,Wage);

% B. Fast interpolation of K''
Kpp_bot                        = sum(vK<=Kp);
Kpp_bot(Kpp_bot<1)             = 1;
Kpp_bot(Kpp_bot>=size(vK,1))   = size(vK,1)-1;
Kpp_up                         = Kpp_bot+1;
Bot_weight                     = (vKNextGuess(Kpp_up)-Kp/(vKNextGuess(Kpp_up)-vKNextGuess(Kpp_bot)));
Kpp                            = Bot_weight * vKNextGuess(Kpp_bot) + (1 - Bot_weight) * vKNextGuess(Kpp_up);

% C. Get tomorrow's elements
[Np,Cp,~,~]             = fnFindN(Kp,Kpp,BCError,ConsImplied,Interest,Wage);

% D. Get the residual
Res                     = EulerError(C,N,K,Kp,Cp,Np,Kpp);
end 