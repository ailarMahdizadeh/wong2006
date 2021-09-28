function [JN,I0,ts,g,JAext,tA,dt,ONOFF,snoise]=CreateModel_2006()
JN.OO = 0.2609; % nA
JN.TT = JN.OO;
JN.OT= 0.0497; % nA
JN.TO= JN.OT;
I0 = 0.3255; % nA
ts = 0.1; % s %tau_NMDA
g = 0.641; %unitless %gamma
JAext = 0.2243E-3; 
tA = 0.002; % s %tau_AMPA
dt = 0.001; %s

ONOFF=1;
snoise=0.02;