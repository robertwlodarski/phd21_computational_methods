# Main functions

# 1. Aiyagari (1994) general equilibrium function

##################################################
##          1. Aiyagari (1994) GE               ##
##################################################

function fnSolveAiyagari1994(p)

        ## 1.   Prepare the setting
        # Unpack parameters
        @unpack α = p

        #       Function settings 
        iWeightOld          = 0.95
        iErrorGE            = 10.0
        iTolGE              = 1e-8
        iTolVFI             = 1e-8
        iTolDist            = 1e-8
        iAccInterval        = 20
        iAccStart           = 30
        iIterNumGE          = 1

        #       Unpack grids 

        #       Preallocation 

end
