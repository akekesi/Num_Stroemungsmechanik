%% Numeriche Stroemungsmechanik
%  Normierte Raeume
%  Numerik von Andreas Malcherek (YouTube)
%  https://www.youtube.com/watch?v=rX4fsL2SXPA&list=PLeJlNT9hA2Pwn8dEA_oJhoD2xEU9iwYMY&index=3
clear all
clc

% Norm
A = [1 2 3;4 5 6;7 8 -9];

% Spaltensummennorm
ssn(A)
norm(A,1)

% Quadratsummennorm
qsn(A)
%sumsqr(A)

% Zeilensimmennorm
zsn(A)
norm(A,inf)

%% Spaltensummennorm
function [x] = ssn(A)
[m,n] = size(A);
spalt = zeros(1,n);
for n = 1:1:n
    for m = 1:1:m
        spalt(1,n) = spalt(1,n)+abs(A(m,n));
    end
end
x = maximum(spalt);
end

%% Quadratsummennorm
function [x] = qsn(A)
[m,n] = size(A);
x = 0;
for m = 1:1:m
    for n = 1:1:n
        x = x+A(m,n)^2;
    end
end
x = sqrt(x);
end

%% Zeilensummennorm
function [x] = zsn(A)
[m,n] = size(A);
zeile = zeros(1,m);
for m = 1:1:m
    for n = 1:1:n
        zeile(1,m) = zeile(1,m)+abs(A(m,n));
    end
end
x = maximum(zeile);
end

%% Maximum Spaltenvektor
function [x] = maximum(A)
x = A(1,1);
for n = 1:1:length(A)
    if x < A(1,n)
        x = A(1,n);
    end
end
end