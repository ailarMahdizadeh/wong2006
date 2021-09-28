function [s0] = PhasePlane_without_stimulus()


thresh = .5;
NN = 5;
%coher = 0;
NFF = 1;
% TCohs = [6.4]'./100;
% coher = repmat(TCohs,NFF,1);

 [JN,I0,ts,g,JAext,tA,dt,ONOFF,snoise]=CreateModel_2006();
function d = dS1(S1, S2, I0, I1, Inoise1)
    d = -S1/ts + (1-S1).*g.*H(Isyn1(S1,S2, I0, I1, Inoise1));
end

function d = dS2(S1, S2, I0, I2, Inoise2)
    d = -S2/ts + (1-S2).*g.*H(Isyn2(S1, S2, I0, I2, Inoise2));
end

function x = Isyn1(S1, S2, I0, I1, Inoise1)
    x = JN.OO*S1 - JN.OT*S2 + I0 + I1 + Inoise1;
end

function x = Isyn2(S1, S2, I0, I2, Inoise2)
    x = JN.TT*S2 - JN.TO*S1 + I0 + I2 + Inoise2;
end

function h = H(x)
    a = 270; b = 108; d = 0.154;
    h = (a*x - b) ./ (1 - exp(-d*(a*x - b)));
    h(h < 0) = 0; 
end

function d = dInoise(Inoise)
    d = 1/tA * -(Inoise + randn(size(Inoise))*sqrt(tA/dt*snoise^2));
end
I1=0;
I2=0;

% function [I1, I2] = stimulus(coh,mu0,ONOFF) 
%     I1 = JAext * mu0 * (1+coh)* ONOFF;
%     I2 = JAext * mu0 * (1-coh)* ONOFF;
% end
 
  % Hz, the baseline external input firing rate.

%Euler integration:
%Euler integration:
function [t, history] = euler(Func, var0, dt, time0, time);
   
    
    t = time0:dt:time;
    history = zeros([size(var0) numel(t)]);
    history(:,:,1) = var0;
    
    for i = 1:numel(t);
        var0 = var0 + Func(var0) .* dt;
        history(:,:,i) = var0;
    end
end

function dx = step(x)
    dx(:,4) = dInoise(x(:,4));
    dx(:,1) = dS1(x(:,1), x(:,2), I0, I1, x(:,3));
    dx(:,2) = dS2(x(:,1), x(:,2), I0, I2, x(:,4));
    dx(:,3) = dInoise(x(:,3));
end
%% 
figure, hold on;
[S2, S1] = meshgrid(linspace(0, 0.8, 25));
%[I1, I2] = stimulus(coher,mu0,NFF);
% axis([0 0.8 0 0.8]); axis square;
% quiver(S1, S2, dS1(S1, S2, I0, I1, 0), dS2(S1, S2, I0, I2, 0));
ylabel(' (S1,Right)');
xlabel(' (S2,Left)');

%% Plot nullclines and fixed points
%[I1, I2] = stimulus(coher,mu0,NFF);
dS = @(x)[ dS1(x(:,1), x(:,2), I0, I1, 0), ...
          dS2(x(:,1), x(:,2), I0, I2, 0)];

s0 = fsolve(dS, [0.03 0.2;0.2 0.03; 0.2 0.2; 0.1 0.7; 0.7 0.1]);
%s0;
 scatter(s0(:,1), s0(:,2), 100);

S = 0:0.01:0.8;
S10 = zeros(length(S),2);
S20 = zeros(length(S),2);
k=0;
for S = 0:.01:0.8
    k=k+1;    

    f1 = @(x) dS1(S, x, I0, I1, 0);
    f2 = @(x) dS2(x, S, I0, I2, 0);
    S20(k,:)= fsolve(f1, 0);
    S10(k,:)= fsolve(f2, 0);
end;
S = 0:0.01:0.8;
plot(S,S10(:,1),'r');
plot(S20(:,1),S,'g');

ylim([0 0.8]); 
xlim([0 0.8]); 
axis square;
legend boxoff;
set(gcf,'Color','w');
set(gca,'Box','off');
set(gca,'FontSize',20);
title('phase plane without stimulus')
end
