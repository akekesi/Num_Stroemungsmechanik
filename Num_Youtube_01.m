%% Numeriche Stroemungsmechanik
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=1kEigsmKlaE&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=1
clear all
clc

%% Wachstum der Bevoelkerung
global c

c = 0.1;
u0 = 0.001;
tanf = 0;
tend = 100;
dt = 2;
t = tanf:dt:tend;
Lt = length(t);

u_ana = u0*exp(c*t);
[u_Eu,u_LW,u_3_,u_RK] = berechnung(u0,Lt,dt);

%% Plot
plot(t,u_ana,'k')
hold on
plot(t,u_Eu,'Color','#0072BD')
plot(t,u_LW,'Color','#D95319')
plot(t,u_3_,'Color','#EDB120')
plot(t,u_RK,'Color','#7E2F8E')
legend({'ana','expl. Euler','Lax-Wendroff','-3-','Runge-Kutta 4'})
grid on
grid minor
hold off

%% Berechnung der Werte
function [u_Eu,u_LW,u_3_,u_RK] = berechnung(u0,Lt,dt)
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

for n = 1:1:Lt-1
    u_Eu(1,n+1) = u_Eu(1,n)+dt*A(u_Eu(1,n));
    u_LW(1,n+1) = u_LW(1,n)+dt*A(u_LW(1,n))+1/2*dt*A(A(u_LW(1,n)));
    u_3_(1,n+1) = u_3_(1,n)+dt*A(u_3_(1,n))+1/2*dt*A(A(u_3_(1,n)))+1/6*dt^2*A(A(A(u_3_(1,n))));
    u_RK(1,n+1) = u_RK(1,n)+dt*A(u_RK(1,n))+1/2*dt*A(A(u_RK(1,n)))+1/6*dt^2*A(A(A(u_RK(1,n))))+1/24*dt^3*A(A(A(A(u_RK(1,n)))));
end
end

%% A-Operator
function [Au] = A(u)
global c
Au = c*u;
end