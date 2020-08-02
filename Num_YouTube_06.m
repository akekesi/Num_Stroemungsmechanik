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
xe = 2*pi;
dx = (xe-xa)/(N-1);
x = xa:dx:xe;

% Zeit
ta = 0;
te = 1;

a = 1;
Ne = 0.4;   % Neumann-Zahl <= 1/2 --> exp. Verf. stabil 
dt = Ne*dx^2/a;

%u(1,:) = rand(1,N);
u(1,:) = sin(x+1/4*pi);

% rb_d Randbed. Dirichlet
% rb_n Randbed. Neumann
f_1 = @(t) 0*t;
f_N = @(t) 1*t;

[u_expl_d] = expl(u,a,dx,ta,te,dt,"rb_d");
[u_impl_d] = impl(u,a,dx,ta,te,dt,"rb_d");
[u_cn_d] = cn(u,a,dx,ta,te,dt,"rb_d");
[u_expl_n] = expl(u,a,dx,ta,te,dt,"rb_n",f_1,f_N);
[u_impl_n] = impl(u,a,dx,ta,te,dt,"rb_n",f_1,f_N);
[u_cn_n] = cn(u,a,dx,ta,te,dt,"rb_n",f_1,f_N);
[u_mat] = matlabdel2(u,a,ta,te,dt);

[steps,~] = size(u_expl_d);

%% Plot
% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_06.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Hilfswerte
ymax = max(u(1,:));
yMax = fix(15*ymax)/10;
ymin = min(u(1,:));
yMin = fix(15*ymin)/10;

% Plot
for m = 1:1:steps
    subplot(2,1,1)
    p11 = plot(x,u_expl_d(1,:),'k');
    hold on
    p_ana = plot(x,x,'Color','#77AC30');
    p12 = plot(x,u_expl_d(m,:),'-x','Color','#0072BD');
    p13 = plot(x,u_impl_d(m,:),'-o','Color','#EDB120');
    p14 = plot(x,u_cn_d(m,:),'-+','Color','#D95319');
    p_mat = plot(x,u_mat(m,:),'--','Color','#4DBEEE');
    p_del2 = plot(x,4*del2(u(1,:),x),':');
    title({'Fourier-Waermeleitung/Diffusion','RB: Dirichlet'},'FontSize',16,'FontWeight','normal')
    legend([p11 p12 p13 p14 p_mat p_del2 p_ana ],{'f(x,t=0)','f(x,t) - expl.','f(x,t) - impl.','f(x,t) - CN','?del2?','2.Abl.','?ana?'},'location','NorthEast')
    xlim([xa,xe])
    ylim([yMin,yMax])
    xlabel('x')
    ylabel('f(x,t)')
    text(xe*0.8,ymax*0.85,['t = ',num2str(round((ta+(m-1)*dt),2)),'s'])
    text(xe*0.8,ymax*0.65,['Ne = ',num2str(Ne)])
    grid on
    grid minor
    drawnow
    hold off
    
    subplot(2,1,2)
    p11 = plot(x,u_expl_n(1,:),'k');
    hold on
    p_ana = plot(x,x,'Color','#77AC30');
    p12 = plot(x,u_expl_n(m,:),'-x','Color','#0072BD');
    p13 = plot(x,u_impl_n(m,:),'-o','Color','#EDB120');
    p14 = plot(x,u_cn_n(m,:),'-+','Color','#D95319');
    p_del2 = plot(x,4*del2(u(1,:),x),':');
    title({'Fourier-Waermeleitung/Diffusion','RB: Neumann'},'FontSize',16,'FontWeight','normal')
    legend([p11 p12 p13 p14 p_del2 p_ana ],{'f(x,t=0)','f(x,t) - expl.','f(x,t) - impl.','f(x,t) - CN','2.Alb.','?ana?'},'location','NorthEast')
    xlim([xa,xe])
    ylim([yMin,yMax])
    xlabel('x')
    ylabel('f(x,t)')
    text(xe*0.8,ymax*0.85,['t = ',num2str(round((ta+(m-1)*dt),2)),'s'])
    text(xe*0.8,ymax*0.65,['Ne = ',num2str(Ne)])
    grid on
    grid minor
    drawnow
    hold off

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

%% Matlab del2
%  https://www.youtube.com/watch?v=snRYb2yG8BQ
%  ??? DE MIT CSINAL A RANBENDINGUNGOKKAL ???
function [u] = matlabdel2(u,a,ta,te,dt)

m = 2;
t = ta+dt;
while t < te
    u(m,:) = u(m-1,:)+a*del2(u(m-1,:));
    m = m+1;
    t = t+dt;
end
end
%% expl
function [u] = expl(u,a,dx,ta,te,dt,RB,funk_1,funk_N)
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
while t < te
    u(m,:) = (M+dt*D)*u(m-1,:)';
    u(m,1) = u(m,1)-dt/dx*funk_1(t-dt);     % ??? -dt ???
    u(m,end) = u(m,end)+dt/dx*funk_N(t-dt); % ??? -dt ???
    m = m+1;
    t = t+dt;
end
end

%% impl
function [u] = impl(u,a,dx,ta,te,dt,RB,funk_1,funk_N)
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
while t < te
    u(m-1,1) = u(m-1,1)-dt/dx*funk_1(t);        % ??? t ???
    u(m-1,end) = u(m-1,end)+dt/dx*funk_N(t);    % ??? t ???
    u(m,:) = (M-dt*D)\u(m-1,:)';
    m = m+1;
    t = t+dt;
end
end

%% Crank-Nicolson
function [u] = cn(u,a,dx,ta,te,dt,RB,funk_1,funk_N)
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
while t < te
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