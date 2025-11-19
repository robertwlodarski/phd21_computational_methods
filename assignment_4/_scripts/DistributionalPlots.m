%% Plotting distributions

% Find maximum index, a = 60
Distributions.pLimitA       = sum(Grids.vGridA2<60);

% Colours
Distributions.vColours      = orderedcolors("gem");

figure(11);
plot(Grids.vGridA2(1:Distributions.pLimitA),cumsum(Results.mMarginalWealth(1,1:Distributions.pLimitA)),'LineWidth',2,'Color',Distributions.vColours(1,:));
hold on;
grid on;
plot(Grids.vGridA2(1:Distributions.pLimitA),cumsum(Results.mMarginalWealth(50,1:Distributions.pLimitA)),'LineWidth',2,'Color',Distributions.vColours(1,:),'LineStyle','--');
plot(Grids.vGridA2(1:Distributions.pLimitA),cumsum(ResultsPE.mMarginalWealth(1,1:Distributions.pLimitA)),'LineWidth',2,'Color',Distributions.vColours(2,:));
plot(Grids.vGridA2(1:Distributions.pLimitA),cumsum(ResultsPE.mMarginalWealth(50,1:Distributions.pLimitA)),'LineWidth',2,'Color',Distributions.vColours(2,:),'LineStyle','--');

