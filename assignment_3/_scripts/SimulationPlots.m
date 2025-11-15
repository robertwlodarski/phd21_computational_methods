% Plot simulations

%% 1. All simulated agents

figure(1);
plot(1:1:10000,mean(Simulations.mWealth,2),'HandleVisibility','off');