%% Numeriche Stroemungsmechanik
%  Transport+Diffusion
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=14
clear all
clc
clf

vid = 0;        % als Video zu speichern: vid = 1

% Konst.
CFL = 0.8;      % bei 0.9 instabil nach ein best. t  
nue = 0.025;

% Gebiet
N = 500;
xa = 0;
x0 = 5;
xe = 20;
dx = (xe-xa)/N;
x = linspace(xa+dx/2,xe-dx/2,N);

% Geschw.
u0 = 4.91;
u = u0*ones(size(x));

% Zeit
ta = 0;
t0 = 1;
te = 2;
dt = 0.01;
dt = zeitschritt(dt,dx,CFL,max(u)); % u = const.
t   = ta:dt:te;

c_ana = @(x,t) 1/sqrt(4*nue*(t0+t))*exp(-(x-x0-u0*t).^2/(4*nue*(t0+t)));
[C_ana] = ana(c_ana,x,t);
c0 = C_ana(1,:);
[C_fvm_upstr] = fvm_upstr(C_ana(1,:),u,nue,dx,N,t,dt,C_ana(1,1),C_ana(1,end));

%% Plot
% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_11.avi');
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
    p11 = plot(x,C_ana(1,:),'Color','[0.7 0.7 0.7]');
    hold on
    p_ana = plot(x,C_ana(m,:),'Color','k');
    p_fvm = plot(x,C_fvm_upstr(m,:),'Color','#0072BD');
    ylim([yMin yMax])
    title({'Transport+Diffusion','1D-FVM-Upstream'},'FontSize',16,'FontWeight','normal')
    legend([p11 p_ana p_fvm],{'f(x,t=0)','f(x,t) - ana.','f(x,t) - FVM-Upstream',},'location','NorthEast')
    text(xe*0.92,ymax*1.00,['CFL = ',num2str(CFL)])
    text(xe*0.92,ymax*0.97,['dt = ',num2str(round(dt,4))])
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

%% Analytisch
function [C] = ana(c_ana,x,t) %,x0,t,t0,u0,nue
C(1,:) = c_ana(x,t(1));
for m = 2:1:length(t)
    C(m,:) = c_ana(x,t(m));
end
end

%% Upstream
function [C] = fvm_upstr(c,u,nue,dx,N,t,dt,RB1,RB2)
C(1,:) = c;

for m = 2:1:length(t)
    for n = 1:1:N
        if n == 1
            ucp = (u(n)+abs(u(n)))/2*C(m-1,n)+(u(n+1)-abs(u(n+1)))/2*C(m-1,n+1);
            fp = ucp-nue/dx*(C(m-1,n+1)-C(m-1,n));
            fn = RB1;
        elseif n == N
            ucn = (u(n-1)+abs(u(n-1)))/2*C(m-1,n-1)+(u(n)-abs(u(n)))/2*C(m-1,n);
            fp = RB2;
            fn = ucn-nue/dx*(C(m-1,n)-C(m-1,n-1));
        else
            ucp = (u(n)+abs(u(n)))/2*C(m-1,n)+(u(n+1)-abs(u(n+1)))/2*C(m-1,n+1);
            ucn = (u(n-1)+abs(u(n-1)))/2*C(m-1,n-1)+(u(n)-abs(u(n)))/2*C(m-1,n);
            fp = ucp-nue/dx*(C(m-1,n+1)-C(m-1,n));
            fn = ucn-nue/dx*(C(m-1,n)-C(m-1,n-1));
        end
        C(m,n) = C(m-1,n)-dt/dx*(fp-fn);
    end
end
end

%% Zeitschritt
function [dt] = zeitschritt(dt,dx,CFL,u)
dt = min(max(CFL*dx/u,10^-14),dt);
end