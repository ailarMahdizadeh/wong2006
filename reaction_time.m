clc;
clear;
close all;
%% Run Model
NumIter = 10;
TCohs = [0 3.2 6.4 12.8 25.6 51.2]'./100;
Coh = repmat(TCohs,NumIter,1);
RW = 1; 
IW = 1;
mu0 =30;
snoise = 0.02;
thresh = 0.4; 
MuLims = linspace(10,40,10);
for i=1:length(MuLims);
    [t, history] = DynamicModel(Coh,RW,IW,MuLims(i),snoise); 
    [ACC(i),RT(i)]=GetBehave(history,thresh);
end

figure, hold on;
% plot(MuLims,RT*1000,'k.-','LineWidth',2,'MarkerSize',20);
% xlabel('Mu0');
% ylabel('RT (ms)');
% set(gcf,'Color','w');
% set(gca,'Box','off');
% set(gca,'FontSize',20);
plot(MuLims,ACC,'k.-','LineWidth',2,'MarkerSize',20);
xlabel('Mu0');
ylabel('ACC');
set(gcf,'Color','w');
set(gca,'Box','off');
set(gca,'FontSize',20);