##          Inner functions

#           Content:
#           1.      Log-normal Tauchen
#           2.      MMV wealth grid


##################################################
##          1. Prepare log-normal Tauchen       ##
##################################################

function        fnTauchenLogNormal(p,N,M)
    #           Log-normal structure
    #           X = log(Z)
    #           Compute the grid for X, and then use Z = exp(X)
    #           P remains the same 

    # 1.        Grid preparation
    #           Unpacking parameters
    @unpack σ, ρ = p

    #           Preparing the x grid
    X_min       = - N * σ / (sqrt(1 - ρ^2))
    X_max       = - X_min
    vGridX      = collect(range(X_min,X_max,length=M))
    Step        = vGridX[2]-vGridX[1]
    vGrid       = exp.(vGridX)

    # 2.        Transition probabilities
    mTransition = zeros(M,M)
    for j in 1:M
        for i in 1:M
            IntPlus             = (vGridX[j] + Step/2  - ρ * vGridX[i]) / σ
            IntMinus            = (vGridX[j] - Step/2  - ρ * vGridX[i]) / σ
            Prob                = cdf(Normal(),IntPlus)-cdf(Normal(),IntMinus)
            mTransition[i,j]    =Prob
        end
    end

    # 3. Print the results
    return vGrid, mTransition
end 

##################################################
##          2. Prepare MMV wealth grid          ##
##################################################

function        fnGridMMV(Min,Max,Points)
    # Finer grid near lower wealth (Maliar, Maliar, and Valli, 2010)
    x           = collect(range(0,0.5,length=Points))
    y           = x.^5 / maximum(x.^5)
    vGrid       = Min .+ (Max - Min) .* y 
    return vGrid
end 