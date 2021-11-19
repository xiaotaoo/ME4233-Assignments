Ni_x=15; Ni_y=5;    
Lx=1; Ly=1;                 %define the length of the domain
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
   
%QR Decomposition
[soln_QR]= QR(A,b);

QR_reshape=reshape(soln_QR,Ni_x-1,Ni_y-1); %reshape into a 2d matrix 
QR_reshape=[0             zeros(1,Ni_y-1) 0
       zeros(Ni_x-1,1) QR_reshape         zeros(Ni_x-1,1)
       0             zeros(1,Ni_y-1) 0 ];  %change into the actual grid with the boundary conditions
   
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
loglog(hori, error_list, 'b');
hold on
loglog(hori, 2*hori.^2, '-r');
legend("Numerical", "RL2 = Ch2")
hold off