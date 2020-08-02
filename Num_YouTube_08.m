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

% Gebiet NxN Knoten
N = 100;
xa = 0;
xe = 1;
h = (xe-xa)/(N-1);
x = xa:h:xe;
ya = xa;
ye = xe;
y = ya:h:ye;

% Zeit
ta = 0;
steps = 100+1;

% Eisen
rho = 7871;
c = 439;
lambda = 80.2;

a = lambda/(rho*c);
Ne = 0.1;   % Neumann-Zahl <= 1/2 --> exp. Verf. stabil 
dt = Ne*h^2/a;

u(:,:,1) = anfangswerte(N);

% rb_d Randbed. Dirichlet
% rb_n Randbed. Neumann
f_1 = @(t) 0*t;
f_N = @(t) 0*t;

u_del2 = matlabdel2(u,steps);
u_expl = expl_Euler(u,a,x,y,h,dt,steps);

%% Plot
% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_08.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Plot
for m = 1:1:steps
    subplot(1,2,1)
    colormap(hot)
    X = [xa xe];
    Y = [ya ye];
    image(X,Y,u_del2(:,:,m))
    colorbar
    title({'Fourier-Waermeleitung/Diffusion','MATLAB del2'},'FontSize',16,'FontWeight','normal')
    text(0.8,0.1,'?dt?')
    text(0.8,0.15,['Steps = ',num2str(m-1)])
    daspect([1 1 1])
    drawnow

    subplot(1,2,2)
    colormap(hot)
    X = [xa xe];
    Y = [ya ye];
    image(X,Y,u_expl(:,:,m))
    colorbar
    title({'Fourier-Waermeleitung/Diffusion','expl. Euler'},'FontSize',16,'FontWeight','normal')
    text(0.8,0.1,['?Ne = ',num2str(Ne),'?'])
    text(0.8,0.15,['dt = ',num2str(round(dt,3))])
    text(0.8,0.2,['Steps = ',num2str(m-1)])
    text(0.8,0.25,['t = ',num2str(round((dt*(m-1)),2))])
    daspect([1 1 1])
    drawnow

    if vid == 1
        frame = getframe(gcf);
        writeVideo(Video,frame);
    end
end

if vid == 1
    close(Video)
end

%% Matlab del2
function [u] = matlabdel2(u,steps)
for s = 2:1:steps
    u(:,:,s) = u(:,:,s-1) + del2(u(:,:,s-1));
end
end

%% expl. Euler
function [u] = expl_Euler(u,a,x,y,h,dt,steps)
nx = length(x);
ny = length(y);
A = laplace(nx,ny);
rb = rb_d(nx,ny,u(:,:,1));
s = 2;
while s <= steps
    u(:,:,s) = u(:,:,s-1)+dt*reshape(a/h^2*(A*reshape(u(:,:,s-1),nx*ny,1)+rb),ny,nx);
    s = s+1;
end
end

%% Laplace
function [A] = laplace(nx,ny)
A = zeros(nx*ny);
for m = 1:1:nx*ny
    A(m,m) = -4;
    if mod(m,ny) == 1
        A(m,m+1) = 1;
        A(m+1,m) = 1;
    elseif mod(m,ny) == 0
        A(m,m-1) = 1;
        A(m-1,m) = 1;
    else
        A(m,m-1) = 1;
        A(m,m+1) = 1;
    end
    if m <= (nx-1)*ny
        A(m,m+ny) = 1;
        A(m+ny,m) = 1;
    end
end
A = sparse(A);
end

%% Randbed.: Dirichlet
function [rb] =rb_d(nx,ny,u)
RB = zeros(ny+2,nx+2);
RB(1,2:end-1) = u(1,:);
RB(end,2:end-1) = u(end,:);
RB(2:end-1,1) = u(:,1);
RB(2:end-1,end) = u(:,end);
rb = zeros(ny*nx,1);
rb(1,1) = RB(1,2)+RB(2,1);
rb(ny,1) = RB(end,2)+RB(end-1,1);
rb((nx-1)*ny+1,1) = RB(1,end-1)+RB(2,end);
rb(end,1) = RB(end,end-1)+RB(end-1,end);
for m = 2:1:(nx-1)*ny
    if m < ny
        rb(m,1) = RB(m+1,1);
        rb((nx-1)*ny+m,1) = RB(m+1,nx+2);
    elseif mod(m,ny) == 1
        rb(m,1) = RB(1,(m-1)/ny+2);
        rb(m+ny-1,1) = RB(ny+2,(m-1)/ny+2);
    end
end
end

%% Anfangswerte
function [u] = anfangswerte(N)
u = ones(N,N)*75;
for m = 1:1:N
    for n = 1:1:N
        if (m-N*0.50)^2+(n-N*0.50)^2 < N
            u(m,n) = 150;
        elseif m > N*0.25 && m < N*0.35 && n > N*0.25 && n < N*0.35
            u(m,n) = 50;
        elseif m > N*0.70 && m < N*0.80 && n > N*0.50 && n < N*0.90
            u(m,n) = 200;
        elseif m > N*0.50 && m < N*0.80 && n > N*0.70 && n < N*0.80
            u(m,n) = 350;
        end
    end
end
end