function [Result]       = fnSoftLog(X)
    if X>0
        Result = log(X);
    elseif X<=0 
        Result = -1e9;
    end
end