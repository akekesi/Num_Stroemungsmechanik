%% Numeriche Stroemungsmechanik
%  FEM - Poisson-GL
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=19

clear all
clc
clf

% field/mesh
xa = -10;
xe = 10;
dx = 0.1;
x = xa:dx:xe;
N = length(x);

% matrices/operators
M = 4*diag(ones(N,1));
D = -2*diag(ones(N,1));
for n = 1:1:N-1
    M(n,n+1) = 1;
    M(n+1,n) = 1;
    D(n,n+1) = 1;
    D(n+1,n) = 1;
end

f = sin(x);
b = M*f'*dx^2/6;
ana = -sin(x)+1/10*(3+sin(10)-2.5)*x+2.5;

% boundary condition - Dirichlet
uL = 2;
uR = 3;
D(1,:) = 0;
D(1,1) = 1;
D(N,:) = 0;
D(N,N) = 1;
b(1) = uL;
b(N) = uR;

% solution
u = D\b;

% plot
plot(x,f)
hold on
plot(x,u)
plot(x,ana,'.k')
legend({'f(x)','u(x)','u_{ana}(x)'},'location','NorthEast')
title({'Poisson-GL-1D: u" = f'},'FontSize',16,'FontWeight','normal')
grid on
