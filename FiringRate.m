clc;
clear;
close all;
%% 
NumIter = 2; %%10
thresh_gate = 0.5;
thresh_fire=15;
% TCohs = [0 3.2 6.4 12.8 25.6 51.2]'./100;
TCohs = [0 51.2]'./100;
Coh = repmat(TCohs,NumIter,1);
C_0 = {'b','r'}; % Cell array of colros.
C=repmat(C_0,1,NumIter);
%% Run Model
[t, history,firing1,firing2] = SimpleModel(Coh);

%% Plot

figure, hold on;
for i = 1:length(Coh)
    
    if squeeze(history(i,1,end))>thresh_gate % Winner (1)
       P{1}= plot(t, squeeze(history(i, 1, :)), 'b');
    else
        plot(t, squeeze(history(i, 1, :)), 'b--');
    end;
    if squeeze(history(i,2,end))>thresh_gate % Loser (2)
        plot(t, squeeze(history(i, 2, :)), 'r--');
    else
       P{2}= plot(t, squeeze(history(i, 2, :)), 'r');
    end;
    
end
plot(t,thresh_gate*ones(size(t)),'k--','LineWidth',1)
legend([P{1} P{2}],{'Pop1','Pop2'});
legend boxoff
ylim([0 0.8]); axis square;
xlabel('Time (s)');
ylabel('S (a.u.)');
set(gcf,'Color','w');
set(gca,'Box','off');
set(gca,'FontSize',20);
%% Coherence Level Effect on Firing Rate of Winner Only
figure, hold on;

for i = 1:length(Coh)
    TmpCoh = find(TCohs==Coh(i));
    Color = [1-TmpCoh/length(TCohs) 0 TmpCoh/length(TCohs)];
    if squeeze(history(i,1,end))>thresh_gate
        plot(t, squeeze(history(i, 1, :)), 'Color',Color);
    else
        plot(t, squeeze(history(i, 1, :)),'--','Color',Color);
    end;
    
    
end
plot(t,thresh_gate*ones(size(t)),'k--','LineWidth',1)
ylim([0 0.8]); axis square;
xlabel('Time (s)');
ylabel('S (a.u.)');
set(gcf,'Color','w');
set(gca,'Box','off');
set(gca,'FontSize',20);
%% Fring rates
figure, hold on;
for i = 1:length(Coh)
    
    if firing1(i,end)>thresh_fire % Winner (1)
       P{1}= plot(t, firing1(i,:), 'color',C{i});
    else
        plot(t, firing1(i,:),'--', 'color',C{i});
    end;
    if firing2(i,end)>thresh_fire % Loser (2)
        plot(t, firing2(i,:), 'color',C{i});
    else
       P{2}= plot(t, firing2(i,:),'--', 'color',C{i});
    end;
    legendInfo{i} = ['Coh = ' num2str(Coh(i))]; 
end
AA=unique(legendInfo);
legend(legendInfo{1:4});

plot(t,thresh_fire*ones(size(t)),'k--','LineWidth',1)
% legend(['c1=', num2str(TCohs(1)),' %'],['c2=', num2str(TCohs(2)),' %'])
% legend([Tcohs(1), Tcohs(2)],{'coh1','coh2'})
% legend([C{1} C{2}],{'c1','c2'});
% legend([P{1} P{2}],{'pop1','pop2'});
legend boxoff
ylim([0 25]); axis square;
xlabel('Time (s)');
ylabel('Firing Rates (Hz)');
set(gcf,'Color','w');

set(gca,'FontSize',20);
%%  Firing Rate of Winner Only
figure, hold on;

for i = 1:length(Coh)
    TmpCoh = find(TCohs==Coh(i));
    Color = [1-TmpCoh/length(TCohs) 0 TmpCoh/length(TCohs)];
    if firing1(i,end)>thresh_fire
        plot(t, firing1(i,:), 'Color',Color);
    else
        plot(t, firing1(i,:),'--','Color',Color);
    end;
    
    
end
plot(t,thresh_fire*ones(size(t)),'k--','LineWidth',1)
ylim([0 50]); axis square;
xlabel('Time (s)');
ylabel('Firing rates (Hz)');

set(gcf,'Color','w');
set(gca,'Box','off');
set(gca,'FontSize',20);
