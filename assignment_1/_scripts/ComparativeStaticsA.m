% Comparative statics (A)

% Original
[w_star, ~]    = fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize);

% Define & discretise the function
Wage              = @(A) arrayfun(@(A_scalar) fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,A_scalar,pa,pr,pb,pMaxIter,pStepSize),A);


% Plots
figure(5); 
fplot(Wage, [0.5 1.5], 'LineWidth', 2.5);
yline(w_star,'--', '$w^\star$','interpreter','latex', 'LabelVerticalAlignment','bottom','LineWidth',2.5,'FontSize',18);
grid on;
xlabel('A');
xlabel('w');
%title('Bisection: Excess supply of labour (aggregate)','fontsize',18);
%legend('Excess supply', 'Location', 'best');
%saveas(gcf,'_figures/excess_supply_bisection.png');
