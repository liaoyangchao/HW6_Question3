%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Homework 6, Question 3
%%%% Student Name: Yangchao Liao
%%%% Student ID.: 1299252
%%%% Department: Civil & Environmental Eng.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
close all;
clc;

%% Initial and boundary conditions
% Define the domain
h = 0.1;
x = 0: h: 2;
y = 0: h: 1;
M = length(x);
N = length(y);
dx = x(2) - x(1);
dy = y(2) - y(1); %new

% Set initial condition for u
u = zeros(M,N);
for i=1:N
    u(1,i) = 0;  % u(x = 0,y) = 0
end
u(M,:) = y'; % u(x = 2,y) = y

u_Point_Jacobi = u;
for i = 2: M-1
    for j = 2: N-1
        u_Point_Jacobi(i,j) = ((1/(dx*dx)) * ( u_Point_Jacobi(i-1,j) + u_Point_Jacobi(i+1,j)) + (1/(dy*dy)) * ( u_Point_Jacobi(i,j-1) + u_Point_Jacobi(i,j+1)))/(2*((1/(dx*dx))+1/(dy*dy)));
    end
end

figure(1)
[X,Y] = meshgrid(y,x);
surf(X,Y,u_Point_Jacobi);hold on;
shading faceted
colormap(cool)

%% Exact solution
x_one = ones(1,11);
n = 1;
u_exact = 0;
if n<Inf
    u_exact = u_exact-4*(1/((n*pi)*(n*pi)*sinh(2*n*pi)))*sinh(n*pi*x)'.*cosh(n*pi*y);
    n = n + 2;
end

u_exact =  (x'/4)*x_one-u_exact;
  
a = x'/4;
b = (x'/4)*x_one;

figure(2)
[X,Y] = meshgrid(y,x);
surf(X,Y,u_exact);hold on;
shading faceted
colormap(cool)