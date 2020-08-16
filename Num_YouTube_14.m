%% Numeriche Stroemungsmechanik
%  FEM - Diffusionsgleichung
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=20

%  PROBLEM: am Anfang gibt es relativ grosse Abweichung mit t=1 (in YouTube auch)
%           --> t=0.1 --> kleinere Abweichung

clear all
clc
clf

vid = 1;        % als Video zu speichern: vid = 1

% field/mesh
xa = 0;
xe = 50;
dx = 0.1;
x = (xa:dx:xe)';
nx = length(x);

% time
t0 = 1;
dt = 0.1;
nt = 100;

% data
nue = 1;
Ne = nue*dt/dx^2;

% matrices/operators
A = (2/3+2*Ne)*diag(ones(nx,1));
B = 4*diag(ones(nx,1));
for n = 1:1:nx-1
    A(n,n+1) = 1/6-Ne;
    A(n+1,n) = 1/6-Ne;
    B(n,n+1) = 1;
    B(n+1,n) = 1;
end

c_ana = @(x,t) 1/sqrt(4*nue*(t+t0))*exp(-(x-25).^2/(4*nue*(t+t0)));
c0 = c_ana(x,0);
c = c0;

% boundary condition - Dirichlet
A(1,:) = 0;
A(1,1) = 1;
A(nx,:) = 0;
A(nx,nx) = 1;

% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_14.avi');
    Video.FrameRate = 13;
    open(Video)
end

% solution
for n = 1:1:nt
    b = B*c/6;
    b(1) = c_ana(x(1),n*dt);
    b(nx) =c_ana(x(end),n*dt);
    c = A\b;

    % plot
    plot(x,c0)
    hold on
    plot(x,c)
    plot(x,c_ana(x,n*dt),'.k')
    legend({'c(x,t = 0)','c(x,t)','u_{ana}(x,t)'},'location','NorthEast')
    title({'Diffusionsgleichung-1D'},'FontSize',16,'FontWeight','normal')
    ylim([0 0.7])
    grid on
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
