function [A,b,r]=construct_grid(Nx,Ny,dx,dy)
%assuming we use the second order central finite difference scheme
%(0.25)(d^2u/dx^2) + (1/9)(d^2u/dy^2) = -(21)(8)(pi^2)cos[(2pi)(2x+3y)]


dx2=dx*dx;              %define dx^2 to be used in the equation
dy2=dy*dy;              %define dy^2 to be used in the equation
dx2dy2=dx2*dy2;         %define dx^2*dy^2

A=zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1));   %we are constructing the internal grid points (square matrix A), so we need to minus 1
b=zeros((Nx-1)*(Ny-1),1);               %we are constructing b from the source term (RHS)of the poisson equation
r=zeros((Nx-1)*(Ny-1),1);


for j=2:Ny-2                            %in this loop, we are iterating through the points not associated with the boundary
    for i=2:Nx-2                        %all the elements in the grid are not in contact with the edge, so they are defined
        po=i+(j-1)*(Nx-1);          % when we loop twice, we are essentially flattening the grid into one column in the x direction then the y direction
        A(po,po)= -2*(9*dy2+4*dx2);    % append the coefficient of (i,j)
        A(po,po+1)=9*dy2;             % append the coefficient of (i+1,j)
        A(po,po-1)=9*dy2;             % append the coefficient of (i-1,j)
        A(po,po-(Nx-1))=4*dx2;        % append the coefficient of (i,j-1)
        A(po,po+(Nx-1))=4*dx2;        % append the coefficient of (i,j+1)
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy))); %this is g(i,j) * (dx^2)(dy^2) at the associated positions
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy))); %obtain a vector r to use for part d
    end
end

% Bottom of the grid
j=1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
        A(po,po-1)=9*dy2;
%         A(po,po-(Nx-1))=4*dx2;      %Since this is the bottom of the grid,
                                    %the row below will be in contact with
                                    %the edge, which must be separately defined
                                    %by the boundary condition u=0
        A(po,po+(Nx-1))=4*dx2; 
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% Top of the grid
j=Ny-1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
        A(po,po-1)=9*dy2;
        A(po,po-(Nx-1))=4*dx2;
%         A(po,po+(Nx-1))=4*dx2; Since this is the top of the grid, the row
                                %below will be in contact with the edge
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% Right side of the grid
i=1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
%         A(po,po-1)=dy2; likewise, on the right edge, there is the
                            % boundary where u=0
        A(po,po-(Nx-1))=4*dx2;
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% Left side of the grid
i=Nx-1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
%         A(po,po+1)=dy2; replace with boundary condition of u=0
        A(po,po-1)=9*dy2;
        A(po,po-(Nx-1))=4*dx2;
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% Bottom left edge of the grid
i=1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
%         A(po,po-1)=dy2;       at the corners, there are two adjacent
%         A(po,po-(Nx-1))=dx2;  nodes that are at the edge(boundary)
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));

% Bottom right edge of the grid  
i=Nx-1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
%         A(po,po+1)=9*dy2; 
        A(po,po-1)=9*dy2;
%         A(po,po-(Nx-1))=4*dx2;
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));

% Top right edge of the grid    
i=Nx-1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
%         A(po,po+1)=9*dy2; 
        A(po,po-1)=9*dy2;
        A(po,po-(Nx-1))=4*dx2;
%         A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));

% Top left edge of the grid        
i=1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
%         A(po,po-1)=9*dy2;
        A(po,po-(Nx-1))=4*dx2;
%         A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        
end