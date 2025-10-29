% Comparative statics (A)

% Original
[w_star, ~]    = fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,pA,pa,pr,pb,pMaxIter,pStepSize);

% Define & discretise the function
Wage              = @(A) arrayfun(@(A_scalar) fnSolvePrices(pTau,pEta,pChi,pBeta,pSigmaz,pAlpha,A_scalar,pa,pr,pb,pMaxIter,pStepSize),A);
AggregateLabSup   = @(w,T) arrayfun(@(w_scalar,T_scalar) fnAggregateLabourSupply(w_scalar,T_scalar,pTau,pEta,pChi,pBeta,pSigmaz,pa,pr,pb,pMaxIter,pStepSize),w,T);
pAmin             = 0.5;
pAmax             = 1.5; 
vSpaceA           = linspace(pAmin, pAmax, 10);
[vWageA, vTA]     = Wage(vSpaceA);
[~,vEmpA,~]       = AggregateLabSup(vWageA,vTA);

% Plot options
Plot5.vPointsX    = [1,1.25];
Plot5.vPointsY    = Wage(Plot5.vPointsX);
Plot5.vPointAnn   = []


% Plot
figure(5); 
plot(vSpaceA,vWageA, 'LineWidth', 2.5,'Color',"b");
hold on;
xlim([pAmin pAmax]);
yline(w_star,'--', '$w^\star$','interpreter','latex', 'LabelVerticalAlignment','bottom','LineWidth',2.5,'FontSize',18);
grid on;
xlabel('A');
ylabel('w');
plot(Plot5.vPointsX,Plot5.vPointsY,'r*');
%title('Bisection: Excess supply of labour (aggregate)','fontsize',18);
%legend('Excess supply', 'Location', 'best');
%saveas(gcf,'_figures/excess_supply_bisection.png');
