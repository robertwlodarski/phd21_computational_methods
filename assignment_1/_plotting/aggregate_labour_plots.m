% Aggregate plots (visualising bisection)

%% 1. Prepare plot values

% Solve the model
[w_star, T_star]    = fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize);

% Define RHS & LHS 
vTTemp              = T_star;
vzTemp              = pAvgz;

% Create RHS and LHS values
Supply              = @(w) arrayfun(@(w_scalar) fnAggregateLabourSupply(w_scalar, vTTemp, pTau, pEta, pChi, pBeta, pSigmaz, pa, pr, pb, pMaxIter, pStepSize), w);
Demand              = @(w) arrayfun(@(w_scalar) fnLabourDemandInverse(pa,w_scalar,pA,pAlpha), w);
ExcessSupply        = @(w) Supply(w) - Demand(w);

% Plot
figure(2); 
fplot(ExcessSupply, [0.3 1], 'LineWidth', 2.5);
xline(w_star,'--', '$w^\star$','interpreter','latex', 'LabelVerticalAlignment','bottom','LineWidth',2.5,'FontSize',18);
xline(w_star-0.05,'-.', '$w_{\text{lower}}$','interpreter','latex', 'LabelVerticalAlignment','bottom','LineWidth',1.5, 'Color','c','FontSize',18);
xline(w_star+0.05,'-.', '$w_{\text{upper}}$','interpreter','latex', 'LabelVerticalAlignment','bottom','LineWidth',1.5, 'Color','c','FontSize',18);
yline(0,'-','Color','r','LineWidth',2.5)
grid on;
xlabel('w');
%title('Bisection: Excess supply of labour (aggregate)','fontsize',18);
legend('Excess supply', 'Location', 'best');
saveas(gcf,'_figures/excess_supply_bisection.png');

