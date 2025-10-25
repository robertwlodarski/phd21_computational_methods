% Plotting working households labour supply curve 


%% 1. Prepare plot values

% Equilibrium values
[w_star, T_star]    = fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize);

% Define RHS & LHS 
vwTemp              = w_star;
vTTemp              = T_star;
vzTemp              = pAvgz;

% Create RHS and LHS values
fnRHS               = @(c) fnIntensiveLabourSupplyRHS(vwTemp,c,vTTemp,vzTemp,pTau,pEta,pChi,pa,pr);
fnLHS               = @(c) (1+pBeta) * c;
fnRHSlow            = @(c) fnIntensiveLabourSupplyRHS(0.8*vwTemp,c,vTTemp,vzTemp,pTau,pEta,pChi,pa,pr);

%% 2. Plot

figure(1); 
fplot({fnRHS,fnRHSlow, fnLHS}, [0.3 2], 'LineWidth', 2.5);
grid on;
legend({'RHS ($w = w^\star$)','RHS ($w = 0.8w^\star$)', 'LHS'},...
    'fontsize',14,'interpreter','latex','Location','best');
% title('Working household Euler equation: Consumption decision','fontsize',18);
% subtitle('Assumed values: $T = 2$, $z = 1$, $\chi = 1$',...
%     'fontsize',14,'interpreter','latex');
xlabel('c');
saveas(gcf,'_figures/consumption_plot.png');
