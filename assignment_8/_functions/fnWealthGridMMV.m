function vGrida      = fnWealthGridMMV(Min, Max, NumberPoints)
    % Finer grid near lower wealth (Maliar, Maliar, and Valli, 2010)
    x               = linspace(0,0.5,NumberPoints)';
    y               = x.^5/max(x.^5);
    vGrida          = Min + (Max - Min) * y;
end