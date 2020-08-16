%% Numeriche Stroemungsmechanik
%  FEM - Lam. Grenzschicht
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=21

%  PROBLEM:

clear all
clc
clf

vid = 1;        % als Video zu speichern: vid = 1

% field/mesh
za = 0;
ze = 1;
dz = 0.01;
z = (za:dz:ze)';
nz = length(z);

% time
dt = 0.01;
nt = 100;

% data
nue = 1;
Ne = nue*dt/dz^2;
gJ = 0.01;

% matrices/operators
A = (2/3+2*Ne)*diag(ones(nz,1));
B = 4*diag(ones(nz,1));
J = z./z;
J(1) = 1/2;
J(nz) = 0;
for n = 1:1:nz-1
    A(n,n+1) = 1/6-Ne;
    A(n+1,n) = 1/6-Ne;
    B(n,n+1) = 1;
    B(n+1,n) = 1;
end

% boundary condition - Dirichlet
A(1,:) = 0;
A(1,1) = 1;

% boundary condition - Neumann
A(nz,nz) = 1/3+Ne;
B(nz,nz) = 2;

% initial condition
v = 0*z;

% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_15.avi');
    Video.FrameRate = 13;
    open(Video)
end

% solution
for n = 1:1:nt
    b = B*v/6+gJ*dt/dz^2*J;
    % boundary condition - Dirichlet
    b(1) = 0;
    v = A\b;

    % plot
    plot(v,z)
    legend({'v_x(z,t)'},'location','NorthEast')
    title({'Lam. Grenzschicht'},'FontSize',16,'FontWeight','normal')
    xlabel('v_x(z,t)')
    ylabel('z')
    xlim([0 50])
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
