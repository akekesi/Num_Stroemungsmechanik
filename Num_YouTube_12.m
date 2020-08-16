%% Numeriche Stroemungsmechanik
%  Burgersgleichung
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=17

%  PROBLEM bei ???

clear all
clc
clf

vid = 0;        % als Video zu speichern: vid = 1

% Field
xa = -2;
xe = 2;
nx = 100;
dx = (xe-xa)/nx;
x = linspace(xa,xe,nx+1);           % edges (x = xa:dx:xe)
xh = linspace(xa+dx/2,xe-dx/2,nx);  % inside

% Time
nt = 100;     % Numer of time steps
ta = 0;
te = 10;

% Dates
g = 9.81;
nue = 0.05;
gamma = 0.1;    % friction coefficient/resistance/Reibungsbeiwert
CFL = 0.8;

% Initial conditions
z_w = 0.5+0.5*exp(-xh.^2/0.2^2);    % height of water
z_b = 0.2+0.0*xh.^2;                % height of ground
u = zeros(size(x));
h = max(z_w-z_b,0);
hu = zeros(size(xh));

% Boundary condition
uL = 0;
uR = 0;

%{
subplot(2,1,1)
plot(xh,z_b,'k')
hold on
plot(xh,z_w,'Color','#0072BD')

subplot(2,1,2)
plot(x,u,'k')
%}

% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_12.avi');
    Video.FrameRate = 13;
    open(Video)
end

t = ta;
for T = 1:1:nt
    % Zeitschritt
    dt = min(0.01,max(CFL/(max(abs(u))+2*nue/dx)/dx,10^(-14)));

    % Bestimmung von h an u-Knoten fuer Durchfluss
    hu(1) = max([0 h(1)]);
    for n = 2:1:nx
        hu(n) = max([0 h(n-1) h(n)]);
    end
    hu(nx+1) = max([0 h(nx)]);
    
    % Gammastar
    for n = 1:1:nx+1
        gammastar(n) = 1+dt*gamma*abs(u(n))./max(hu(n),0.001);
    end

    % ??? Fu = fvm_upst(u,u,dx,dt,0,uL,uR,CFL)/gammastar(n); %(n); ???
    Fu = fvm_upst(u,u,dx,dt,0,uL,uR,CFL);
    for n = 1:1:nx+1
        Fu(n) = Fu(n)/gammastar(n);
    end

    % freesurf()
    for n = 1:1:nx
        if n == 1
            lo(n) = 0;
            up(n) = -g*dt^2/dx^2*hu(n+1);
            di(n) = 1-lo(n)-up(n);
            b(n) = z_w(n)-dt/dx*(hu(n+1)*Fu(n+1)-hu(n)*uL);
        elseif n == nx
            lo(n) = -g*dt^2/dx^2*hu(n);
            up(n) = 0;
            di(n) = 1-lo(n)-up(n);
            b(n) = z_w(n)-dt/dx*(-hu(n+1)*uR+hu(n)*Fu(n));
        else
            lo(n) = -g*dt^2/dx^2/gammastar(n)*hu(n);     % ????? -+ ????? video--><--pdf
            up(n) = -g*dt^2/dx^2/gammastar(n+1)*hu(n+1); % ????? -+ ????? video--><--pdf
            di(n) = 1-lo(n)-up(n);
            b(n) = z_w(n)-dt/dx*(hu(n+1)*Fu(n+1)-hu(n)*Fu(n));
        end
    end
    zss = thomas(lo,di,up,b);
    z_w = max(zss',z_b);
    h = max(z_w-z_b,0);

    % LSG der Imp-GL
    for n = 1:1:nx+1
        if n == 1
            u(n) = uL;
        elseif n == nx+1
            u(n) = uR;
        else
            u(n) = Fu(n)-g*dt/dx*(z_w(n)-z_w(n-1))/gammastar(n);
        end
    end
    if t+dt > te
        dt = te-t;
    end
    if t >= te
        break
    end

    subplot(2,1,1)
    plot(xh,z_b,'k')
    hold on
    plot(xh,z_w,'Color','#0072BD')
    ylim([-1 1])
    grid on
    grid minor
    hold off

    subplot(2,1,2)
    plot(x,u,'k')
    ylim([-1 1])
    grid minor
    hold off

    drawnow
    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

function [c] = fvm_upst(c,u,dx,dt,nue,RB1,RB2,CFL)
nc = length(c);
tau = 0;    % initial local time
dtau = min(max(CFL*dx/max(max(abs(u)),10^(-14))),dt);
cnew = zeros(size(c));
for k = 1:1:1000
    tau = tau+dtau;
    if tau >= dt
        break
    end
    if tau+dtau > dt
        dtau = dt-dtau; % letzter Schritt muss passen
    end
    for n = 1:1:nc
        if n == 1
            ucp = (u(n+1)+abs(u(n+1)))/2*c(n)+(u(n+1)-abs(u(n+1)))/2*c(n+1);
            fp = ucp-nue/dx*(c(n+1)-c(n));
            fn = RB1;
        elseif n == nc
            ucn = (u(n)+abs(u(n)))/2*c(n-1)+(u(n)-abs(u(n)))/2*c(n);
            fp = RB2;
            fn = ucn-nue/dx*(c(n)-c(n-1));
        else
            ucp = (u(n+1)+abs(u(n+1)))/2*c(n)+(u(n+1)-abs(u(n+1)))/2*c(n+1);
            ucn = (u(n)+abs(u(n)))/2*c(n-1)+(u(n)-abs(u(n)))/2*c(n);
            fp = ucp-nue/dx*(c(n+1)-c(n));
            fn = ucn-nue/dx*(c(n)-c(n-1));
        end
        cnew(n) = c(n)-dtau/dx*(fp-fn);
    end
end
c = cnew;
end

function [x] = thomas(al,ac,au,b) % book s73 Numerische Methoden der Stroemungsmechanik
% al: lower diagonal, a(1) is not used
% ac: diagonal
% au: upper diagonal
% b:  right hand side
N = length(b);
au(1) = au(1)/ac(1);
b(1) = b(1)/ac(1);
for n = 2:1:N
    gamma = 1/(ac(n)-au(n-1)*al(n));
    au(n) = au(n)*gamma;
    b(n) = (b(n)-al(n)*b(n-1))*gamma;
end
x = zeros(N,1);
x(N) = b(N);
for n = N-1:-1:1
    x(n) = b(n)-au(n)*x(n+1);
end
end
%{
function [] = freesurf() % ?????????????
for n = 1:1:nx
    if n == 1
        lo(n) = 0;
        up(n) = -g*dt^2/dx^2*hu(n+1);
        di(n) = 1-lo(n)-up(n);
        b(i) = z_w(n)-dt/dx*(hu(n+1)*FU(n+1)-hu(n)*uL);
    elseif n == nx
        lo(n) = -g*dt^2/dx^2*hu(n);
        up(n) = 0;
        di(n) = 1-lo(n)-up(n);
        b(i) = z_w(n)-dt/dx*(-hu(n+1)*uR+hu(n)*Fu(n));
    else
        lo(n) = g*dt^2/dx^2*gammastar(n)*hu(n);
        up(n) = g*dt^2/dx^2*gammastar(n+1)*hu(n+1);
        di(n) = 1-lo(n)-up(n);
        b(i) = z_w(n)-dt/dx*(hu(n+1)*Fu(n+1)-hu(n)*Fu(n));
    end
end
zss = thomas(lo,di,up,b);
z_w = max(zss',z_b);
h = max(z_w-z_b,0);
end
%}