clf%% Numeriche Stroemungsmechanik
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=rX4fsL2SXPA&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=5
clear all
clc

dx = 0.05;
xa = 0;
xe = 1;
x = xa:dx:xe;

u = 2.5;          % Geschw.
nue = 0.1;      % kin. Visk.
Pe = dx*u/nue;  % Pecletzahl
STABILITAET_ZD = 'Stabilitaet (Zentralendiff.): wenn A diagonaldominant ist, d.h. Pe < 2';
STABILITAET_VD = 'Stabilitaet (Vorderendiff.): ??? unabhängig von Pe ??? --> NALAM NEM, nalam Pe<=1-nel stabil';
%% Berechnung

% numerisch
[A_zd,b_zd] = operator_zd(x,Pe);
c_num_zd = A_zd\b_zd;

[A_vd,b_vd] = operator_vd(x,Pe);
c_num_vd = A_vd\b_vd;

% analytisch
c_ana = @(x) (exp(u/nue*x)-1)/(exp(u/nue)-1);

%% Plot

plot(x,c_ana(x),'LineWidth',1.5,'Color','#0072BD')
hold on
plot(x,c_num_zd,'-x','LineWidth',1.5,'Color','#D95319')
plot(x,c_num_vd,'-+','LineWidth',1,'Color','#EDB120')
legend('ana','num (Zentralendiff)','num (Vorderendiff)')
text(0.1,0.9,['Pe = ',num2str(Pe)])
text(0.1,0.85,STABILITAET_ZD)
text(0.1,0.8,STABILITAET_VD)
grid on
grid minor
hold off

%% A - Verfahrensoperator mit Zentralendifferenz
function [A,b] = operator_zd(x,Pe)
N = length(x);
A = zeros(N,N);
b = zeros(N,1);
for m = 1:1:N
    for n = 1:1:N
        if n == m
            A(m,n) = -2;
        elseif n == m-1
            A(m,n) = 1+Pe/2;
        elseif n == m+1
            A(m,n) = 1-Pe/2;
        end
    end
end
% Randbedingung
% c(x=0) = 0;
A(1,:) = 0;
A(1,1) = 1;
% c(x=1) = 1;
A(N,:) = 0;
A(N,N) = 1;
b(n,1) = 1;
end

%% A - Verfahrensoperator mit Vorderendifferenz
function [A,b] = operator_vd(x,Pe)
N = length(x);
A = zeros(N,N);
b = zeros(N,1);
for m = 1:1:N
    for n = 1:1:N
        if n == m
            A(m,n) = 2-Pe;
        elseif n == m-1
            A(m,n) = -1;
        elseif n == m+1
            A(m,n) = Pe-1;
        end
    end
end
% Randbedingung
% c(x=0) = 0;
A(1,:) = 0;
A(1,1) = 1;
% c(x=1) = 1;
A(N,:) = 0;
A(N,N) = 1;
b(n,1) = 1;
end