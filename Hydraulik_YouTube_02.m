%% Hydraulik - YouTube
%  Reservoir
%  https://www.youtube.com/watch?v=Fz8zPAPJLA4
clear all
clc
clf

vid = 1;    % als Video zu speichern: vid = 1

global d0
global d1
global d2
global d3
global h0
global h3
global te
global n_ueber
global h_ueber
global t_ueber

d0 = 5;
d1 = 0.5;
d2 = 0.3;
d3 = 0.3;
A0 = d0^2*pi()/4;
A1 = d1^2*pi()/4;
A2 = d2^2*pi()/4;
A3 = d3^2*pi()/4;
h0 = 4;
h3 = 3;
w1 = 10;
w2 = 1;
w3 = 5;

h = zeros(1)+2;
t = zeros(1);
te = 30;
dt = 0.1;
n = 2;

n_ueber = 1;
h_ueber = [];
t_ueber = [];


% Video Open
if vid == 1
    Video = VideoWriter('Hydraulik_YouTube_02.avi');
    Video.FrameRate = 13;
    open(Video)
end

% Berechnung/Calculatio
while t < te
    if h(n-1) < h3
        h(n) = h(n-1)+(A1*w1-A2*w2)/A0*dt;
    elseif h(n-1) >= h3 && h(n-1) < h0
        h(n) = h(n-1)+(A1*w1-A2*w2-A3*w3)/A0*dt;
    elseif h(n-1) >= h0
        h(n) = h0;
    else
        h(n) = 0;
    end
    t(n) = t(n-1)+dt;
    plot_funk(h,t,vid,Video)
    n = n+1;
end

% Video Close
if vid == 1
    close(Video)
end

%% Plot
function [] = plot_funk(h,t,vid,Video)
global d0
global d1
global d2
global d3
global h0
global h3
global te
global n_ueber
global h_ueber
global t_ueber

L = d0*0.1;
subplot(1,2,1)
plot([-d0/2-L d0/2+L],[0 0],'k','LineWidth',3)
hold on
plot([-d0/2 -d0/2],[d1 h0],'k','LineWidth',3)
plot([-d0/2-L -d0/2],[d1 d1],'k','LineWidth',3)
plot([d0/2 d0/2],[d2 h3],'k','LineWidth',3)
plot([d0/2 d0/2+L],[d2 d2],'k','LineWidth',3)
plot([d0/2 d0/2],[h3+d3 h0],'k','LineWidth',3)
plot([d0/2 d0/2+L],[h3 h3],'k','LineWidth',3)
plot([d0/2 d0/2+L],[h3+d3 h3+d3],'k','LineWidth',3)
if h(end) >= h0
    Col = '#D95319';
else
    Col = '#0072BD';
end
plot([-d0/2 d0/2],[h(end) h(end)],'Color',Col,'LineWidth',5)
ylim([0 h0*1.2])
title({'Animation'},'FontSize',16,'FontWeight','normal')
hold off

subplot(1,2,2)
p = plot(t,h,'k');
hold on
if length(h) > 1 && h(end) >= h3 && h(end-1) < h3
    h_ueber(n_ueber) = h(end);
    t_ueber(n_ueber) = t(end);
    n_ueber = n_ueber+1;
end
for m = 1:1:length(h_ueber)
    plot([t_ueber(m) t_ueber(m)],[0 h_ueber(m)],'--k')
    plot(t_ueber(m),h_ueber(m),'k','Marker','.','MarkerSize',15)    
end
xlim([0 te])
ylim([0 h0*1.2])
title({'h - t'},'FontSize',16,'FontWeight','normal')
legend(p,'h(t)','location','NorthEast')
xlabel('t [s]')
ylabel('h(t) [m]')
grid on
grid minor
hold off

drawnow

% Video Save
if vid == 1
    frame = getframe(gcf);
    writeVideo(Video,frame);
end
end