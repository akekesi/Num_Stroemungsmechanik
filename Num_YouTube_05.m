%% Numeriche Stroemungsmechanik
%  Advektions-GL / Lagrange-Verfahren
%  Numerik von Andreas Malcherek (YouTube) ...22/9
%  https://www.youtube.com/watch?v=n_aNqKZPSO8&list=PLeJlNT9hA2PwyyZKrQMcxMauZcHZBi5Rm
clear all
clc
clf

vid = 1;        % als Video zu speichern: vid = 1

% Gebiet
N = 100;
xa = 0;
xe = 10;
dx = (xe-xa)/(N-1);
x = xa:dx:xe;

% Zeit
ta = 0;
te = 5;

a = 1;
CFL = 0.9;
dt = CFL*dx/a;

u_1(1,:) = anfangsfunk_1(x);  % Anfangswerte 1
u_2(1,:) = anfangsfunk_2(x);  % Anfangswerte 2

u_ana_1 = ana(@anfangsfunk_1,x,ta,te,dt,a);
u_ana_2 = ana(@anfangsfunk_2,x,ta,te,dt,a);

[u_11] = upwind(u_1,ta,te,dt,dx,a);
[u_12] = ftcs(u_1,ta,te,dt,dx,a);
[u_13] = lax_fr(u_1,ta,te,dt,dx,a);
[u_14] = lax_we(u_1,ta,te,dt,dx,a);
[u_15] = impl(u_1,ta,te,dt,dx,a);
[u_16] = cn(u_1,ta,te,dt,dx,a);
[u_17] = lagrange(u_1,ta,te,dt,x,dx,a);
[m,n] = size(u_11);
U_1(:,:,1) = u_11;
U_1(:,:,2) = u_12;
U_1(:,:,3) = u_13;
U_1(:,:,4) = u_14;
U_1(:,:,5) = u_15;
U_1(:,:,6) = u_16;
U_1(:,:,7) = u_17;

[u_21] = upwind(u_2,ta,te,dt,dx,a);
[u_22] = ftcs(u_2,ta,te,dt,dx,a);
[u_23] = lax_fr(u_2,ta,te,dt,dx,a);
[u_24] = lax_we(u_2,ta,te,dt,dx,a);
[u_25] = impl(u_2,ta,te,dt,dx,a);
[u_26] = cn(u_2,ta,te,dt,dx,a);
[u_27] = lagrange(u_2,ta,te,dt,x,dx,a);
U_2(:,:,1) = u_21;
U_2(:,:,2) = u_22;
U_2(:,:,3) = u_23;
U_2(:,:,4) = u_24;
U_2(:,:,5) = u_25;
U_2(:,:,6) = u_26;
U_2(:,:,7) = u_27;

%% Plot
% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_05.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
% Hilfswerte
ymax = max(u_ana_1(1,:));
ymin = min(u_ana_1(1,:));
yMax = fix(ymax*15)/10;
if yMax == 0
    yMax = 1;
end
yMin = fix(ymin*15)/10;
if yMin == 0
    yMin = -1;
end
for m = 1:1:m
    subplot(2,1,1)
    p_ana = plot(x,u_ana_1(m,:),'LineWidth',1.5,'Color','k');
    hold on
    p11 = plot(x,U_1(m,:,1),'-o','Color','#77AC30');
    p12 = plot(x,U_1(m,:,2),'Color','#EDB120');
    p13 = plot(x,U_1(m,:,3),'Color','#D95319');
    p14 = plot(x,U_1(m,:,4),'Marker','.','MarkerSize',15,'Color','#0072BD');
    p15 = plot(x,U_1(m,:,5),'-x','Color','#7E2F8E');
    p16 = plot(x,U_1(m,:,6),'-+','Color','#A2142F');
    p17 = plot(x,U_1(m,:,7),'-*','Color','r');
    ylim([yMin yMax])
    title({'Advektions-GL','glatte Anfangsbedingung'},'FontSize',16,'FontWeight','normal')
    legend([p_ana p11 p12 p13 p14 p15 p16 p17],{'ana','Upwind','FTCS','Lax-Friedrich','Lax-Wendroff','impl','Crank-Nicolson','Lagrange'},'location','NorthEast')
    xlabel('x')
    ylabel('f(x)')
    text(xe*0.83,ymax*0.9,['CFL = ',num2str(CFL)])
    text(xe*0.83,ymax*0.75,['t = ',num2str(round((ta+(m-1)*dt),2)),'s'])
    text(xe*0.83,ymax*0.6,['N = ',num2str(N)])
    grid on
    grid minor
    hold off

    subplot(2,1,2)
    p_ana = plot(x,u_ana_2(m,:),'LineWidth',1.5,'Color','k');
    hold on
    p11 = plot(x,U_2(m,:,1),'-o','Color','#77AC30');
    p12 = plot(x,U_2(m,:,2),'Color','#EDB120');
    p13 = plot(x,U_2(m,:,3),'Color','#D95319');
    p14 = plot(x,U_2(m,:,4),'Marker','.','MarkerSize',15,'Color','#0072BD');
    p15 = plot(x,U_2(m,:,5),'-x','Color','#7E2F8E');
    p16 = plot(x,U_2(m,:,6),'-+','Color','#A2142F');
    p17 = plot(x,U_2(m,:,7),'-*','Color','r');
    ylim([yMin yMax])
    title({'Advektions-GL','nichtglatte Anfangsbedingung'},'FontSize',16,'FontWeight','normal')
    legend([p_ana p11 p12 p13 p14 p15 p16 p17],{'ana','Upwind','FTCS','Lax-Friedrich','Lax-Wendroff','impl','Crank-Nicolson','Lagrange'},'location','NorthEast')
    xlabel('x')
    ylabel('f(x)')
    text(xe*0.83,ymax*0.9,['CFL = ',num2str(CFL)])
    text(xe*0.83,ymax*0.75,['t = ',num2str(round((ta+(m-1)*dt),2)),'s'])
    text(xe*0.83,ymax*0.6,['N = ',num2str(N)])
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

%% analytische Loesung
function [u] = ana(anfangsfunk,x,ta,te,dt,a)
u = zeros(1,length(x));
m = 1;
t = ta;
while t < te
    u(m,:) = anfangsfunk(x-a*t);
    m = m+1;
    t = t+dt;
end
end

%% FTCS
function [u] = ftcs(u,ta,te,dt,dx,a)
c = a*dt/dx;

m = 2;
t = ta+dt;
while t < te
    u(m,1) = 0;     % Randbed. 1
    u(m,end) = 0;   % Randbed. 2
    for n = 2:1:length(u(1,:))-1
        u(m,n) = u(m-1,n)-c/2*(u(m-1,n+1)-u(m-1,n-1));
    end
    m = m+1;
    t = t+dt;
end
end

%% Lax-Friedrich
function [u] = lax_fr(u,ta,te,dt,dx,a)
c = a*dt/dx;

m = 2;
t = ta+dt;
while t < te
    u(m,1) = 0;     % Randbed. 1
    u(m,end) = 0;   % Randbed. 2
    for n = 2:1:length(u(1,:))-1
        u(m,n) = 1/2*(u(m-1,n+1)+u(m-1,n-1)-c/2*(u(m-1,n+1)-u(m-1,n-1)));
    end
    m = m+1;
    t = t+dt;
end
end

%% Lax-Wendroff
function [u] = lax_we(u,ta,te,dt,dx,a)
c = a*dt/dx;

m = 2;
t = ta+dt;
while t < te
    u(m,1) = 0;     % Randbed. 1
    u(m,end) = 0;   % Randbed. 2
    for n = 2:1:length(u(1,:))-1
        u(m,n) = u(m-1,n)-c/2*(u(m-1,n+1)-u(m-1,n-1))+c^2/2*(u(m-1,n+1)-2*u(m-1,n)+u(m-1,n-1));
    end
    m = m+1;
    t = t+dt;
end
end

%% Upwind
function [u] = upwind(u,ta,te,dt,dx,a)
ap = (a+abs(a))/2;
an = (a-abs(a))/2;
cp = ap*dt/dx;
cn = an*dt/dx;

m = 2;
t = ta+dt;
while t < te
    u(m,1) = 0;     % Randbed. 1
    u(m,end) = 0;   % Randbed. 2
    for n = 2:1:length(u(1,:))-1
        u(m,n) = u(m-1,n)-cp*(u(m-1,n)-u(m-1,n-1))-cn*(u(m-1,n+1)-u(m-1,n));
    end
    m = m+1;
    t = t+dt;
end
end

%% Implizit
function [u] = impl(u,ta,te,dt,dx,a)
m = 2;
t = ta+dt;
while t < te
    [A,u0] = impl_op(u(m-1,:),dt,dx,a);
    u(m,:) = A\u0;
    m = m+1;
    t = t+dt;    
end
end

%% Implizit-Operator
function [A,b] = impl_op(u,dt,dx,a)
c = a*dt/dx;
N = length(u);
A = zeros(N,N);
b = u';
for m = 1:1:N
    for n = 1:1:N
        if n == m
            A(m,n) = 1;
        end
        if n == m-1
            A(m,n) = -c/2;
        end
        if n == m+1
            A(m,n) = c/2;
        end
    end
end

% Randbedingung
% u(x=0,t) = 0;
A(1,:) = 0;
A(1,1) = 1;
b(1,1) = 0;
% u(x=1,t) = 0;
A(N,:) = 0;
A(N,N) = 1;
b(n,1) = 0;
end

%% Crank-Nicolson
function [u] = cn(u,ta,te,dt,dx,a)
m = 2;
t = ta+dt;
while t < te
    [A,u0] = cn_op(u(m-1,:),dt,dx,a);
    u(m,:) = A\u0;
    m = m+1;
    t = t+dt;    
end
end

%% Crank-Nicolson-Operator
function [A_n1,b] = cn_op(u,dt,dx,a)
theta = 1/2;
c = a*dt/dx;
N = length(u);
A_n = zeros(N,N);
A_n1 = zeros(N,N);
for m = 1:1:N
    for n = 1:1:N
        if n == m
            A_n(m,n) = 1;
            A_n1(m,n) = 1;
        end
        if n == m-1
            A_n(m,n) = (1-theta)*c/2;
            A_n1(m,n) = -theta*c/2;
        end
        if n == m+1
            A_n(m,n) = (theta-1)*c/2;
            A_n1(m,n) = theta*c/2;
        end
    end
end
b = A_n*u';

% Randbedingung
% u(x=0,t) = 0;
A_n1(1,:) = 0;
A_n1(1,1) = 1;
b(1,1) = 0;
% u(x=1,t) = 0;
A_n1(N,:) = 0;
A_n1(N,N) = 1;
b(n,1) = 0;
end

%% Lagrnage
function [u] = lagrange(u,ta,te,dt,x,dx,a)
m = 2;
t = ta+dt;
while t < te
    u(m,1) = 0;     % Randbed. 1
    u(m,end) = 0;   % Randbed. 2
    for n = 2:1:length(u(1,:))-1
        xb = x(n)-a*dt;
        if xb <= x(1)
            u(m,n) = u(m,1);
        elseif xb >= x(end)
            u(m,n) = u(m,end);
        else
            b = (xb-x(1))/dx;
            p = ceil(b);
            d = b-p+1;
            u(m,n) = (1-d)*u(m-1,p)+d*u(m-1,p+1);
        end
    end
    m = m+1;
    t = t+dt;
end
end

%% Anfangsfunktion 1
function [u] = anfangsfunk_1(x)
f = @(x) exp(-5/2*(x-2).^2);
u = f(x);
u(1,1) = 0;        % Randbed. 1
u(1,end) = 0;      % Randbed. 2
end

%% Anfangsfunktion 2
function [u] = anfangsfunk_2(x)
u = zeros(1,length(x));
for n = 1:1:length(x)
    if 1 <= x(1,n) && x(1,n) <= 3
        u(1,n) = 1;
    end
    u(1,1) = 0;     % Randbed. 1
    u(1,end) = 0;   % Randbed. 2
end
end