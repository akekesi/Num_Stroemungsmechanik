%% Numeriche Stroemungsmechanik
%  FEM - Turb. Grenzschicht
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=jIIHfpSVRyY&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=22

%  PROBLEM:

clear all
clc
clf

vid = 1;        % als Video zu speichern: vid = 1

% field/mesh
za = 0;
ze = 1;
dz = 0.02;
z = (za:dz:ze)';
nz = length(z);

% time
dt = 100;
nt = 100;

% data
gJ = 0.2e-6;
ustar = sqrt(gJ*ze);
nuet = 1e-6+ustar*ze*z.*(1-z/ze);
%plot(nuet,z)
%return
Ne = nuet*dt/dz^2;

% matrices/operators
A = 2/3*diag(ones(nz,1));
B = diag(ones(nz,1));
J = z./z;
J(1) = 1/2;
J(nz) = 1/2;
B(1,1) = Ne(1)+1/2*Ne(2);
for n = 2:1:nz-1
    B(n,n) = 1/2*(Ne(n-1)+2*Ne(n)+Ne(n+1));
    B(n,n+1) = -1/2*(Ne(n)+Ne(n+1));
    B(n+1,n) = -1/2*(Ne(n+1)+Ne(n));
end
B(1,2) = -1/2*(Ne(1)+Ne(2));
B(2,1) = -1/2*(Ne(2)+Ne(1));
B(nz,nz) = 1/2*Ne(nz-1)+1/2*Ne(nz);

for n = 1:1:nz-1
    A(n,n+1) = 1/6;
    A(n+1,n) = 1/6;
end
A(nz,nz) = 1/3;

C = A+B;

% boundary condition - Dirichlet
C(1,:) = 0;
C(1,1) = 1;

% initial condition
v = 0*z;

% Video
if vid == 1
    Video = VideoWriter('Num_YouTube_16.avi');
    Video.FrameRate = 13;
    open(Video)
end

% solution
for n = 1:1:nt
    b = A*v+gJ*dt/dz^2*J;
    % boundary condition - Dirichlet
    b(1) = 0;
    v = C\b;

    % plot
    plot(v,z)
    legend({'v_x(z,t)'},'location','NorthEast')
    title({'Turb. Grenzschicht'},'FontSize',16,'FontWeight','normal')
    xlabel('v_x(z,t)')
    ylabel('z')
    xlim([0 10])
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
