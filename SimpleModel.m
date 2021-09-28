function [t, history,firing1,firing2] = SimpleModel(Coh)
mu0=30;
function d = dS1(S1, S2, I0, I1, Inoise1)
    d = -S1/ts + (1-S1).*g.*H(Isyn1(S1,S2, I0, I1, Inoise1));
end

function d = dS2(S1, S2, I0, I2, Inoise2)
    d = -S2/ts + (1-S2).*g.*H(Isyn2(S1, S2, I0, I2, Inoise2));
end

[JN,I0,ts,g,JAext,tA,dt,ONOFF,snoise]=CreateModel_2006();

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



function [I1, I2] = stimulus(coh,mu0,ONOFF) 
    I1 = JAext * mu0 * (1+coh)* ONOFF;
    I2 = JAext * mu0 * (1-coh)* ONOFF;
end




%Euler integration:
function [t, history] = euler(Func, var0, dt, time0, time)    
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

%--------------------------------------------------------------------------

%% Itretion
[I1, I2] = stimulus(Coh,mu0,1);
[t, history] = euler(@step, zeros([length(Coh),4]), dt, 0, 2);

firing1=zeros(numel(Coh),numel(t));
firing2=zeros(numel(Coh),numel(t));


for ii=1:numel(Coh)
firing1(ii,:)=H(Isyn1(history(ii, 1, :),history(ii, 2, :), I0, I1(ii), history(ii, 3, :)));

firing2(ii,:)=H(Isyn2(history(ii, 1, :),history(ii, 2, :), I0, I2(ii), history(ii, 4, :)));
end
end