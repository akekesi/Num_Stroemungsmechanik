%% Numeriche Stroemungsmechanik
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=Ej9Xc3J2GCI&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=2
clear all
clc

%% Wachstum der Bevoelkerung mit logistische DGL
%  Zeitschritt novelessel romlik minden
global gamma
global tau
global theta

umax = 10;
gamma = 0.1;
tau = gamma/umax;
theta = 0.5;
u0 = 0.001;
tanf = 0;
tend = 1000;
dt = 15;
t = tanf:dt:tend;
Lt = length(t);

u_ana = gamma*(tau+(gamma/u0-tau)*exp(-gamma*t)).^(-1);
[u_Eu,u_LW,u_3_,u_RK,u_im,u_CN] = berechnung(u0,Lt,dt);

%% Plot
plot(t,u_ana,'k')
hold on
plot(t,u_Eu,'Color','#0072BD')
plot(t,u_LW,'Color','#D95319')
plot(t,u_3_,'Color','#EDB120')
plot(t,u_RK,'Color','#7E2F8E')
plot(t,u_im,'*','Color','#77AC30')
plot(t,u_CN,'o','Color','#77AC30')
ylim([0 umax*1.2])
legend({'ana','expl. Euler','Lax-Wendroff','-3-','Runge-Kutta 4','impl. Verf.','Crank-Nicolson'})
grid on
grid minor
hold off

%% Berechnung der Werte
function [u_Eu,u_LW,u_3_,u_RK,u_im,u_CN] = berechnung(u0,Lt,dt)
global gamma
global tau
global theta

% expl. Euler
u_Eu = zeros(1,Lt);
u_Eu(1,1) = u0;

% Lax-Wendroff
u_LW = zeros(1,Lt);
u_LW(1,1) = u0;

% 3. ???
u_3_ = zeros(1,Lt);
u_3_(1,1) = u0;

% Runge-Kutta-4
u_RK = zeros(1,Lt);
u_RK(1,1) = u0;

% implizites Verf.
u_im = zeros(1,Lt);
u_im(1,1) = u0;

% Crak-Nicolson.
u_CN = zeros(1,Lt);
u_CN(1,1) = u0;

for n = 1:1:Lt-1
    u_Eu(1,n+1) = u_Eu(1,n)+dt*A(u_Eu(1,n));
    u_LW(1,n+1) = u_LW(1,n)+dt*A(u_LW(1,n))+1/2*dt*A(A(u_LW(1,n)));
    u_3_(1,n+1) = u_3_(1,n)+dt*A(u_3_(1,n))+1/2*dt*A(A(u_3_(1,n)))+1/6*dt^2*A(A(A(u_3_(1,n))));
    u_RK(1,n+1) = u_RK(1,n)+dt*A(u_RK(1,n))+1/2*dt*A(A(u_RK(1,n)))+1/6*dt^2*A(A(A(u_RK(1,n))))+1/24*dt^3*A(A(A(A(u_RK(1,n)))));
    u_im(1,n+1) = -(1-gamma*dt)/(2*tau*dt)+sqrt(((1-gamma*dt)/(2*tau*dt))^2+u_im(1,n)/(tau*dt));
    u_CN(1,n+1) = (-(1-dt*theta*gamma)+sqrt((1-dt*theta*gamma)^2-4*dt*theta*tau*(-1)*(dt*(1-theta)*gamma*u_CN(1,n)-dt*(1-theta)*tau*u_CN(1,n)^2+u_CN(1,n))))/(2*dt*theta*tau);
end
end

%% A-Operator
function [Au] = A(u)
global gamma
global tau
Au = gamma*u-tau*u^2;
end