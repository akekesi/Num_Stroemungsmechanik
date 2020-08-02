%% Hydraulik - YouTube
%  Geschw. im Rohr/Velocity in pipe/tube
%  https://www.youtube.com/watch?v=N80UG9YuLrE
clear all
clc
clf

vid = 0;    % als Video zu speichern: vid = 1

D = 2;
A = D^2*pi()/4;
U = D*pi();
L = 1000;
rho = 1000;
p1 = 105000;
p2 = 100000;
z1 = 0;
z2 = 0;
v0 = 0;
g = 9.81;
lambda = 0.02;

% Zeit/Time
ta = 0;
te = 1000;
dt = 10;
t = ta:dt:te;

% dv Funktion/Function
dv = @(t,v) (p1-p2)/L/rho-g*(z2-z1)/L-lambda/8*U/A*v*abs(v);

%% Berechnung mit RK4/Calculation with RK4
V_rk4 = zeros(1,length(t));
V_rk4(1) = v0;
for n = 2:1:length(t)
    V_rk4(n) = rk4(dv,0,V_rk4(n-1),dt);
end

[P] = visu(D,L,V_rk4,t);

%% Berechnung mit ODE45/Calculation with ODE45
[T,V_ode45] = ode45(dv,[ta te],v0);

%% Plot

% Video
if vid == 1
    Video = VideoWriter('Hydraulik_YouTube_01.avi');
    Video.FrameRate = 13;
    open(Video)
end

for n = 1:1:length(t)
    subplot(2,1,1)
    p_rk4 = plot(t,V_rk4,'Color','k');
    hold on
    p_ode45 = plot(T,V_ode45,'--','Color','#D95319');
    plot(t(n),V_rk4(n),'Marker','.','MarkerSize',15,'Color','k')
    xlabel('t [s]')
    ylabel('v(t) [m/s]')
    legend([p_rk4 p_ode45],{'RK4','ODE45'},'location','NorthEast')
    title({'Geschw. im Rohr'},'FontSize',16,'FontWeight','normal')
    grid on
    hold off

    subplot(2,1,2)

    % Rohr/Pipe/Tube
    plot([0 L],[0 0],'Color',[0 0.4470 0.7410],'LineWidth',120)
    hold on
    plot([0 0],[-D/2 D/2],'k',[L L],[-D/2 D/2],'k',[0 L],[0 0],'-.k')
    plot([0 L],[D/2 D/2],'k','LineWidth',3)
    plot([0 L],[-D/2 -D/2],'k','LineWidth',3)

    % Geschw./Velocity
    vMax = 50;
    vmax = max(V_rk4);
    plot([L/2-vMax/2 L/2+vMax/2],[1.3*D/2 1.3*D/2],'LineWidth',5,'Color',[0.7 0.7 0.7])
    plot(L/2+vMax/2,1.3*D/2,'>','LineWidth',5,'Color',[0.7 0.7 0.7])
    plot([L/2-vMax/2 L/2-vMax/2+V_rk4(n)/vmax*vMax],[1.3*D/2 1.3*D/2],'LineWidth',5,'Color',[0 0.4470 0.7410])
    plot(L/2-vMax/2+V_rk4(n)/vmax*vMax,1.3*D/2,'>','LineWidth',5,'Color',[0 0.4470 0.7410])

    % Pfeile/Arrows
    plot(P(1,:,n),P(2,:,n),'>w')

    xlim([0 L])
    xticks(0:L/10:L)
    ylim([-D D])
    yticks([0 D/2])
    xlabel('L [m]')
    ylabel('R [m]')
    text(L/2-vMax/2,D/2*1.55,['v = ',num2str(round(V_rk4(n),2)),' m/s'])
    text(10,D/2*1.3,['p_1 = ',num2str(round(p1)),' Pa'])
    text(10,D/2*1.7,['z_1 = ',num2str(round(z1)),' m'])
    text(930,D/2*1.3,['p_2 = ',num2str(round(p2)),' Pa'])
    text(930,D/2*1.7,['z_2 = ',num2str(round(z2)),' m'])
    title({'Animation'},'FontSize',16,'FontWeight','normal')
    grid on
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

%% Runge-Kutta-4
function [v_neu] = rk4(funk,t,v,dt)
k1 = funk(t,v);
k2 = funk(t,v+dt/2*k1);
k3 = funk(t,v+dt/2*k2);
k4 = funk(t,v+dt*k3);
v_neu = v+dt*(k1/6+k2/3+k3/3+k4/6);
end

%% Visualisation
function [P] = visu(D,L,V,T)
N = 5;
P = zeros(2,N*N,length(T));
x = L/(N+1)/(N+1);
y = D/2/N;

% Grundpfeil/Basic Arrow
for n = 1:1:N
    P(1,n,1) = n*x;
    P(2,n,1) = D/2-n*y;
    if n ~= N
        P(1,N+n,1) = n*x;
        P(2,N+n,1) = -P(2,n);
    end
end

% Schritte des Grundpfeiles/Steps of Basic Arrow
for n = 2:1:length(T)
    P(1,:,n) = P(1,:,n-1)+V(n-1)*(T(n)-T(n-1));
    P(2,:,n) = P(2,:,n-1);
end

% Pattern/Pattern
for n = 1:1:N
    P(1,n*(2*N-1)+1:(n+1)*(2*N-1),:) = P(1,1:(2*N-1),:)+n*L/(N+1);
    P(2,n*(2*N-1)+1:(n+1)*(2*N-1),:) = P(2,1:(2*N-1),:);
end

% Wiederholung/Replay
for m = 1:1:(2*N-1)*(N+1)
    for n = 1:1:length(T)
        while P(1,m,n) > L
            P(1,m,n) = P(1,m,n)-L;
        end
    end
end
end