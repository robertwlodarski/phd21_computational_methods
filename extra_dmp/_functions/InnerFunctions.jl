# Inner functions

# 1. Plain DMP equilibrium residual

##################################################
##          1. Plain DMP eq. residual           ##
##################################################

function        fnPlainDMPResidual(θ,A,p)
    
    # 1.        Unpack parameters
    @unpack μ,ξ,η,κ,b,β,λ = p

    # 2.        Calculate the implied values 
    q           = μ * θ^(-ξ)
    w           = η*(A + κ *θ) + (1 - η)*b
    J           = (A - w) / (1 - β*(1 - λ))

    # 3.        Get residual 
    return      κ / q - J * β
    
end 