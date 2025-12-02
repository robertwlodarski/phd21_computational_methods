function [a_next, c, n]       = fnFastOptimisation(Wealth,Wage, Interest, LabShock, Minimum,ExpectedValue,Parameters,Grids)
    % Initialisation
    vGridA1             = Grids.vGridA1;
    pBeta               = Parameters.pBeta;
    pMu                 = Parameters.pMu;
    pEta                = Parameters.pEta;
    pChi                = Parameters.pChi;

    % Maximisation function
    Options             = optimset('Display','off','TolFun',1e-12,'TolX',1e-12,'MaxIter',20000,'MaxFunEvals',20000);
    [a_next]            = fminbnd(@fnResidualCalc,Minimum,vGridA1(end),Options);
    A                   = (1 + Interest) * Wealth - a_next - pMu / 2 * (a_next^2 / Wealth - 2*a_next + Wealth);
    B                   = (Wage * LabShock)^2 / pEta;
    c                   = (sqrt(A^2+ 4 * B)+A) / 2;
    n                   = Wage * LabShock / (c * pEta);

    % Nested
    function [OptValue] = fnResidualCalc(Guess)
        % Setting
        a_prime                         = Guess;
        a_bot                           = sum(vGridA1<a_prime);
        a_bot(a_bot <1)                 = 1;
        a_bot(a_bot >= size(vGridA1,1)) = size(vGridA1,1)-1;
        a_up                            = a_bot+1;

        % Weights 
        Weight_low                      = (vGridA1(a_up) - a_prime) / (vGridA1(a_up) - vGridA1(a_bot));
        Weight_low(Weight_low<0)        = 0;
        Weight_low(Weight_low>1)        = 1;

        % Calculations
        AvgExpectedValue                = Weight_low * ExpectedValue(a_bot)+ (1 - Weight_low)*ExpectedValue(a_up);

        % Use the quadratic trick [change unless sigma=chi=1]
        AAA                             = (1 + Interest) * Wealth - a_prime - pMu / 2 * (a_prime^2 / Wealth - 2*a_prime + Wealth); 
        BBB                             = (Wage * LabShock)^2 / pEta;
        Consumption                     = (sqrt(AAA^2+ 4 * BBB)+AAA) / 2;
        Labour                          = Wage * LabShock / (Consumption * pEta);

        % Current period value
        Value                           = log(Consumption)  - pEta / (1 + 1 / pChi) * Labour^(1 + 1 / pChi) + pBeta * AvgExpectedValue;
        OptValue                        = -Value;
    end

end