%% Numeriche Stroemungsmechanik
%  Waermeleitung
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=10
%  https://www.youtube.com/watch?v=snRYb2yG8BQ --> del2 ??? DE MIT CSINAL A RANBENDINGUNGOKKAL ???
clear all
clc
clf

% Fourier-Waermeleitung/Diffusion

vid = 0;        % als Video zu speichern: vid = 1

% Gebiet
N = 100;
xa = 0;
xe = 1;
dx = (xe-xa)/(N-1);
x = xa:dx:xe;

% Zeit
ta = 0;
steps = 10+1;

% Eisen
rho = 7871;
c = 439;
lambda = 80.2;

a = lambda/(rho*c);
Ne = 0.4;   % Neumann-Zahl <= 1/2 --> exp. Verf. stabil 
dt = Ne*dx^2/a;

u(1,:) = anfangswerte(x);

% rb_d Randbed. Dirichlet
% rb_n Randbed. Neumann
f_1 = @(t) 0*t;
f_N = @(t) 0*t;

[u_d] = expl(u,a,dx,ta,steps,dt,"rb_d");
%[u] = impl(u,a,dx,ta,steps,dt,"rb_d");
%[u] = cn(u,a,dx,ta,steps,dt,"rb_d");
[u_n] = expl(u,a,dx,ta,steps,dt,"rb_n",f_1,f_N);
%[u] = impl(u,a,dx,ta,steps,dt,"rb_n",f_1,f_N);
%[u] = cn(u,a,dx,ta,steps,dt,"rb_n",f_1,f_N);

%% Plot
% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_07.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Hilfswerte
[M,~] = size(u_d);

ymax = max(u(1,:));
yMax = fix(ymax*12)/10;
ymin = min(u(1,:));
yMin = fix(ymin*10-abs(ymin*2))/10;

% Plot
for m = 1:1:M
    subplot(2,2,1)
    plot(x,u_d(1,:),'k')
    hold on
    plot(x,u_d(m,:))
    xlim([xa xe])
    ylim([yMin yMax])
    title({'Fourier-Waermeleitung/Diffusion','RB: Dirichlet'},'FontSize',16,'FontWeight','normal')
    legend({'f(x,t=0)','f(x,t) - expl.'},'location','NorthEast')
    text(xe*0.8,ymax*0.95,['Steps = ',num2str(m-1)])
    text(xe*0.8,ymax*0.87,['Ne = ',num2str(Ne)])
    xlabel('Laenge')
    ylabel('Temperatur')
    grid on
    grid minor
    drawnow
    hold off

    subplot(2,2,3)
    colormap(hot)
    Xa = [0 xe];
    Xe = [0 xa];
    image(Xa,Xe,u_d(m,:))
    xlim([xa xe])
    yticks([])
    colorbar
    drawnow
    
    subplot(2,2,2)
    plot(x,u_n(1,:),'k')
    hold on
    plot(x,u_n(m,:))
    xlim([xa xe])
    ylim([yMin yMax])
    title({'Fourier-Waermeleitung/Diffusion','RB: Neumann'},'FontSize',16,'FontWeight','normal')
    legend({'f(x,t=0)','f(x,t) - expl.'},'location','NorthEast')
    text(xe*0.8,ymax*0.95,['Steps = ',num2str(m-1)])
    text(xe*0.8,ymax*0.87,['Ne = ',num2str(Ne)])
    xlabel('Laenge')
    ylabel('Temperatur')
    grid on
    grid minor
    drawnow
    hold off

    subplot(2,2,4)
    colormap(hot)
    Xa = [0 xe];
    Xe = [0 xa];
    image(Xa,Xe,u_n(m,:))
    xlim([xa xe])
    yticks([])
    colorbar
    drawnow

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

%% Anfangswerte
function [u] = anfangswerte(x)
N = length(x);
u = ones(1,N)*50;
for n = 1:1:N
    if n >= 3/5*N && n <= 4/5*N
        u(1,n) = 150;
    elseif n >= 1/5*N && n <= 2/5*N
        u(1,n) = 250;
    end
end
end

%% expl
function [u] = expl(u,a,dx,ta,steps,dt,RB,funk_1,funk_N)
if nargin < 8
  funk_1 = @(t) 0;
  funk_N = @(t) 0;
end

N = length(u);

if RB == "rb_d"
    [M,D] = operator_d(a,dx,N);
elseif RB == "rb_n"
    [M,D] = operator_n(a,dx,N);
end

m = 2;
t = ta+dt;
while m <= steps
    u(m,:) = (M+dt*D)*u(m-1,:)';
    u(m,1) = u(m,1)-dt/dx*funk_1(t-dt);     % ??? -dt ???
    u(m,end) = u(m,end)+dt/dx*funk_N(t-dt); % ??? -dt ???
    m = m+1;
    t = t+dt;
end
end

%% impl
function [u] = impl(u,a,dx,ta,steps,dt,RB,funk_1,funk_N)
if nargin < 8
  funk_1 = @(t) 0;
  funk_N = @(t) 0;
end

N = length(u);

if RB == "rb_d"
    [M,D] = operator_d(a,dx,N);
elseif RB == "rb_n"
    [M,D] = operator_n(a,dx,N);
end

m = 2;
t = ta+dt;
while m <= steps
    u(m-1,1) = u(m-1,1)-dt/dx*funk_1(t);        % ??? t ???
    u(m-1,end) = u(m-1,end)+dt/dx*funk_N(t);    % ??? t ???
    u(m,:) = (M-dt*D)\u(m-1,:)';
    m = m+1;
    t = t+dt;
end
end

%% Crank-Nicolson
function [u] = cn(u,a,dx,ta,steps,dt,RB,funk_1,funk_N)
if nargin < 8
  funk_1 = @(t) 0;
  funk_N = @(t) 0;
end

N = length(u);

if RB == "rb_d"
    [M,D] = operator_d(a,dx,N);
elseif RB == "rb_n"
    [M,D] = operator_n(a,dx,N);
end

m = 2;
t = ta+dt;
theta = 1/2;
while m <= steps
    RB_t = zeros(N,1);
    RB_t(1,1) = dt/dx*funk_1(t);            % ??? dt-t, t ???
    RB_t(end,1) = -dt/dx*funk_N(t);         % ??? dt-t, t ???
    RB_tdt = zeros(N,1);
    RB_tdt(1,1) = dt/dx*funk_1(t-dt);       % ??? dt-t, t ???
    RB_tdt(end,1) = -dt/dx*funk_N(t-dt);    % ??? dt-t, t ???    
    u(m,:) = (M-dt*theta*D)\((M+dt*(1-theta)*D)*u(m-1,:)'-theta*RB_t-(1-theta)*RB_tdt);
    m = m+1;
    t = t+dt;
end
end

%% Operator
%  Randbedingung: Dirichlet
function [M,D] = operator_d(a,dx,N)
M = diag(ones(N,1));
D = diag(ones(N,1))*(-2*a/dx^2);
d1 = ones(N-1,1)*(a/dx^2);
D1 = diag(d1,1);
D2 = diag(d1,-1);
D = D+D1+D2;
D(1,:) = zeros(N,1);    % RB 1
D(end,:) = zeros(N,1);  % RB 2
end

%  Randbedingung: Dirichlet
function [M,D] = operator_n(a,dx,N)
M = diag(ones(N,1));
D = diag(ones(N,1))*(-2*a/dx^2);
d1 = ones(N-1,1)*(a/dx^2);
D1 = diag(d1,1);
D2 = diag(d1,-1);
D = D+D1+D2;
D(1,1) = -a/dx^2;       % RB 1
D(end,end) = -a/dx^2;   % RB 2
end