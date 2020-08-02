%% Numeriche Stroemungsmechanik
%  Transport+Diffusion
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=11
clear all
clc
clf

vid = 1;        % als Video zu speichern: vid = 1

% Gebiet
xa = 0;
x0 = 10;
xe = 100;

% Zeit
ta = 1;
te = 25;

u = 3;
nue = 1;
Pe = 0.9;   % Pe <= 2   Peclet-Zahl
Ne = 0.4;   % Ne <= 1/2 Neumann-Zahl
Cr = 0.9-2*Ne;    % Cr+2Ne <= 1 Courant-Zahl
theta = 1/2;

dx = Pe*nue/u;
dt_Ne = Ne*dx^2/nue;
dt_Cr = Cr*dx/u;
dt = min(dt_Ne,dt_Cr);
x = xa:dx:xe;
t = ta:dt:te;

c_ana = @(x,t) 1/sqrt(4*nue*(t+ta)).*exp(-(x-x0-u*t).^2/(4*nue*(t+ta)));
[C_ana] = ana(c_ana,x,t);                       % Analytisch
[C_ftcs] = ftcs(C_ana(1,:),x,dx,t,dt,u,nue);    % FTCS - NEGATIVE NUMERISCHE DIFFUSION
[C_upstr] = upstr(C_ana(1,:),x,dx,t,dt,u,nue);  % Upstream - ZU VIEL NUMERISCHE DIFFUSION
[C_cn] = cn(C_ana(1,:),x,dx,t,dt,u,nue,theta);
[M,D1,D2,D3,D4] = operator(length(x));

%% Plot
% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_09.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Hilfswert
ymax = max(C_ana(1,:));
yMax = fix(ymax*12)/10;
ymin = min(C_ana(1,:));
yMin = fix(ymin*10-abs(ymin*2))/10;

% Plot
for m = 1:1:length(t)
    p11 = plot(x,C_ana(1,:),'Color',[0.7 0.7 0.7],'LineWidth',2);
    hold on
    p_ana = plot(x,C_ana(m,:),'k','LineWidth',2);
    p_ftcs = plot(x,C_ftcs(m,:),'-x');
    p_upstr = plot(x,C_upstr(m,:),'-+');
    p_cn = plot(x,C_cn(m,:),'-o');
    ylim([yMin yMax])
    title('Transport+Diffusion','FontSize',16,'FontWeight','normal')
    legend([p11 p_ana p_ftcs p_upstr p_cn],{'f(x,t=0)','f(x,t) - ana.','f(x,t) - FTCS','f(x,t) - Upstream','f(x,t) - Crank-Nicolson'},'location','NorthEast')
    text(xe*0.92,ymax*0.95,['Ne = ',num2str(Ne)])
    text(xe*0.92,ymax*0.92,['Pe = ',num2str(Pe)])
    text(xe*0.92,ymax*0.89,['Cr = ',num2str(Cr)])
    text(xe*0.92,ymax*0.86,['\Deltat = ',num2str(dt)])
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

%% Analytitsche
function [C] = ana(funk,x,t)
C(1,:) = funk(x,t(1));
for m = 2:1:length(t)
    C(m,:) = funk(x,t(m));
    C(m,1) = C(1,1);       % Randbed 1
    C(m,end) = C(1,end);   % Randbed 2
end
end

%% FTCS - NEGATIVE NUMERISCHE DIFFUSION
function [C] = ftcs(c0,x,dx,t,dt,u,nue)
[M,D1,D2,~,~] = operator(length(x));
C(1,:) = c0;
for m = 2:1:length(t)
    C(m,:) = (M+dt*nue/dx^2*D1-dt*u/(2*dx)*D2)*C(m-1,:)';
    C(m,1) = c0(1,1);       % Randbed 1
    C(m,end) = c0(1,end);   % Randbed 2
end
end

%% Upstream - ZU VIEL NUMERISCHE DIFFUSION
function [C] = upstr(c0,x,dx,t,dt,u,nue)
[M,D1,~,D3,D4] = operator(length(x));

C(1,:) = c0;
for m = 2:1:length(t)
    C(m,:) = (M+dt*nue/dx^2*D1-dt*(u+abs(u))/(2*dx)*D3-dt*(u-abs(u))/(2*dx)*D4)*C(m-1,:)';
    C(m,1) = c0(1,1);       % Randbed 1
    C(m,end) = c0(1,end);   % Randbed 2
end
end

%% Crank-Nicolson
function [C] = cn(c0,x,dx,t,dt,u,nue,theta)
[M,D1,D2,~,~] = operator(length(x));

C(1,:) = c0;
for m = 2:1:length(t)
    C(m,:) = (M+dt*u/(2*dx)*theta*D2-dt*nue/dx^2*theta*D1)\((M-dt*u/(2*dx)*(1-theta)*D2+dt*nue/dx^2*(1-theta)*D1)*C(m-1,:)');
    C(m,1) = c0(1,1);       % Randbed 1
    C(m,end) = c0(1,end);   % Randbed 2
end
end

%% Operator
function [M,D1,D2,D3,D4] = operator(N)
M = diag(ones(N,1));
d1 = ones(N-1,1);
D1 = diag(ones(N,1))*-2;
D1 = D1+diag(d1,1)+diag(d1,-1);
D2 = zeros(N)+diag(d1,1)-diag(d1,-1);
D3 = diag(ones(N,1))-diag(d1,-1);
D4 = -diag(ones(N,1))+diag(d1,-1);
end