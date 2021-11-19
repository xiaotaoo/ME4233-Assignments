%ME4233 Part 2 Assignment 1
%Finite difference method to solve the poisson equation
%(0.25)(d^2u/dx^2) + (1/9)(d^2u/dy^2) = (21)(8)(pi)cos[(2pi)(2x+3y)]

clear
close all                   %clear and close all previously loaded figures within the command

debug=0;                    %set to 0 to run the actual grid and 1 to debug the code

if debug==0
    Ni_x=60; Ni_y=50;           %define how many grid points in the x and y direction
end
if debug==1
    Ni_x=5; Ni_y=5;           %define how many grid points in the x and y direction
end
Lx=1; Ly=1;                 %define the length of the domain

%Part A: Discretizing the poisson equation using the central difference
x=linspace(0,Lx,Ni_x+1);    %returns a vector of 61 equally spaced elements, 61 elements = 60 grid spaces
dx= x(2)-x(1);              %find the grid space length by subtracting the first element from the second element
y=linspace(0,Ly,Ni_y+1);    %returns a vector of 51 equally spaced elemetns, 51 elements = 50 grid spaces
dy= y(2)-x(1);              %find the grid space length by subtracting the first element from the second element


[A,b,r] = construct_grid(Ni_x,Ni_y,dx,dy);    %call the function named construct grid with our defined arguments

N= size(A,1);               %obtain the integer N which corresponds to the number of rows and columns of A
%Part B: LU and QR decomposition

%LU Decomposition
[soln_LU]= LU(A,b);

LU_reshape=reshape(soln_LU,Ni_x-1,Ni_y-1); %reshape into a 2d matrix 
LU_reshape=[0             zeros(1,Ni_y-1) 0
       zeros(Ni_x-1,1) LU_reshape         zeros(Ni_x-1,1)
       0             zeros(1,Ni_y-1) 0 ];  %change into the actual grid with the boundary conditions
figure(1)
contourf(x,y,LU_reshape')                        %plot the contour plot for the grid
figure(2)
surf(x,y,LU_reshape')


%QR Decomposition
[soln_QR]= QR(A,b);

QR_reshape=reshape(soln_QR,Ni_x-1,Ni_y-1); %reshape into a 2d matrix 
QR_reshape=[0             zeros(1,Ni_y-1) 0
       zeros(Ni_x-1,1) QR_reshape         zeros(Ni_x-1,1)
       0             zeros(1,Ni_y-1) 0 ];  %change into the actual grid with the boundary conditions
figure(3)
contourf(x,y,QR_reshape')                        %plot the contour plot for the grid
figure(4)
surf(x,y,QR_reshape')

%Part C
%Jacobi Method
u0=randi([0,100],N,1); %initialize a random column matrix of size N x 1 for the GS method
[residual_j, iter_j]=jacobi(A,b, soln_QR,u0);

%GS Method
[residual_g, iter_g]=gs(A,b, soln_QR,u0);

%SOR method
L_S=tril(A,-1);     %take the lower triangular matrix of A excluding the diagonal
U_S=triu(A,1);      %take the upper triangular matrix of A excluding the diagonal
iter_s=0;           %initialize the iteration for S
residual_s=[];      %initialize the residual array
u0_s=randi([0,100],N,1);    %randomize a first guess for the SOR method
[residual_s, iter_s]=sor(A,b,soln_QR,u0,1.2);
figure(5)
semilogy(0:iter_j, residual_j, '-*b', 0:iter_g, residual_g, '-*r', 0:iter_s, residual_s, '-*g');
legend('Jacobi Method', 'Gauss-Seidel Method', 'SOR Method, w=1.2')

for w= [0.3 0.5 0.8 1.0 1.2]
    [residual_s, iter_s]=sor(A,b, soln_QR,u0,w);
    disp(['Iteration:' num2str(w)])

    figure(6)
    semilogy(0:iter_s, residual_s)
    hold on
    u0_s=randi([0,100],N,1); 
    iter_s=0;
    residual_s=[];
end
hold off        %stop holding the graphs
legend('w=0.3','w=0.5','w=0.8','w=1.0','w=1.2')

%part d
r1 = randi([0,1000],N,1);   %initializing different random vectors to test
r2 = randi([0,100],N,1);
r3= zeros(N,1);
for r_list = [r r1 r2 r3] 
    [residual_j, iter_j]=jacobi(A,b, soln_QR,r_list);   %try the different random vectors for jacobi, taking the QR soln
    figure (7);
    semilogy(0:iter_j, residual_j);     %plotting syntax
    hold on
end
hold off
legend('r=21*cos(2pi(2x+3y)', 'random 1', 'random 2', 'zero')

for r_list = [r r1 r2 r3] 
    [residual_g, iter_g]=gs(A,b, soln_QR,r_list);  %try the different random vectors for GS, taking the QR soln
    figure (8);
    semilogy(0:iter_g, residual_g);     %plotting syntax
    hold on
end
hold off
legend('r=21*cos(2pi(2x+3y)', 'random 1', 'random 2', 'zero')

for r_list = [r r1 r2 r3] 
    [residual_s, iter_s]=sor(A,b, soln_QR,r_list,w==1.2); %try the different random vectors for SOR, w=1.2
    figure (9);
    semilogy(0:iter_s, residual_s);     %plotting syntax
    hold on
end
hold off
legend('r=21*cos(2pi(2x+3y)', 'random 1', 'random 2', 'zero')

%part e 
error_list = [];
hori = [];  %create a list for the horizontal axis
for Ny =[11 21 31 41 51]
    Nx = Ny + 10;
    left=zeros((Nx-1),1);
    right=zeros((Nx-1),1);
    bottom=zeros(1,Ny-1);
    top=zeros(1,Ny-1);
    analytic = zeros(Nx+1, Ny+1);

    x=linspace(0,Lx,Nx+1);    %returns a vector of 61 equally spaced elements, 61 elements = 60 grid spaces
    dx= x(2)-x(1);              %find the grid space length by subtracting the first element from the second element
    y=linspace(0,Ly,Ny+1);    %returns a vector of 51 equally spaced elemetns, 51 elements = 50 grid spaces
    dy= y(2)-x(1);              %find the grid space length by subtracting the first element from the second element
    [A1,b1,ana] = change_bc(Nx,Ny,dx,dy); %subtract the boundary condition where applicable from the b term
    [soln]= LU(A1,b1);

    error_list= [error_list sqrt(sum(( soln - ana ).^2)/(Ny*Nx))]; % RMS
    hori = [hori sqrt(dx*dy)];
end

figure(10);
loglog(hori, error_list, '-*b');
hold on
loglog(hori, 2*hori.^2, '-r');
legend("Numerical", "RL2 = Ch2")
hold off


