function [t, history] = DynamicModel(Coh,RW,IW,mu0,snoise)

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
    a = 270; b = 108; d = 0.1540;
    h = (a*x - b) ./ (1 - exp(-d*(a*x - b)));
    h(h < 0) = 0; 
end

function d = dInoise(Inoise)
    d = 1/tA * -(Inoise + randn(size(Inoise))*sqrt(tA/dt*snoise^2));
end


tA = 0.002; % s %tau_AMPA
dt = 0.001; %s

function [I1, I2] = stimulus(coh,mu0,ONOFF) 
    I1 = JAext * mu0 * (1+coh)* ONOFF;
    I2 = JAext * mu0 * (1-coh)* ONOFF;
end


JAext = 0.2243E-3; 

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
[t, history] = euler(@step, zeros([length(Coh),4]), dt, 0, 3);
% 
% eq1=-S1/tS + (1-S1).*g.*H(Isyn1(S1,S2, I0, I1, Inoise1));
% eq2 = -S2/tS + (1-S2).*g.*H(Isyn2(S1, S2, I0, I2, Inoise2));
% J=jacobian([eq1,eq2],[S1,S2]);
% x=x_eq;
% J(1,1)= - 3*x^2 + 6*x - 3/50;
% J(2,1)=-10*x;
% disp(J)
% %eigen value
% matrix=(landa*J)-eye(2);
% disp(matrix)
% determinan=det(matrix);
% disp(determinan)
% landa_solve=vpasolve(determinan==0, landa);
% landa=landa_solve(2,1);

end