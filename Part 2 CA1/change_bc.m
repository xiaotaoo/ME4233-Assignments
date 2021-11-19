function [A,b,r]=change_bc(Nx,Ny,dx,dy)
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
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy))); 
    end
end

% Bottom of the grid
j=1;
for i=2:Nx-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
        A(po,po-1)=9*dy2;
%         A(po,po-(Nx-1))=4*dx2*21*cos(4*pi*(0+dx*i));      %this will be shifted to the RHS
        A(po,po+(Nx-1))=4*dx2; 
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-4*dx2*21*cos(4*pi*(0+i*dx));
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
%         A(po,po+(Nx-1))=4*dx2*21*cos(4*pi*(0+dx*i)); Shifted to the RHS

        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-4*dx2*21*cos(4*pi*(0+dx*i));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% Left side of the grid
i=1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
%         A(po,po-1)=9*dy2*21*cos(6*pi*(0+dy*j)); Shifted to the RHS
        A(po,po-(Nx-1))=4*dx2;
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-9*dy2*21*cos(6*pi*(0+dy*j));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% right side of the grid
i=Nx-1;
for j=2:Ny-2
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
%         A(po,po+1)=9*dy2*21*cos(6*pi*(0+dy*j)); replace with boundary condition of u=0
        A(po,po-1)=9*dy2;
        A(po,po-(Nx-1))=4*dx2;
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-9*dy2*21*cos(6*pi*(0+dy*j));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
end

% Bottom left edge of the grid
i=1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
%         A(po,po-1)=9*dy2*21*cos(6*pi*(0+dy*j));       at the corners, there are two adjacent
%         A(po,po-(Nx-1))=4*dx2*21*cos(4*pi*(0+dx*i));  nodes that are at the edge(boundary)
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-4*dx2*21*cos(4*pi*(0+dx*i))-9*dy2*21*cos(6*pi*(0+dy*j));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));

% Bottom right edge of the grid  
i=Nx-1;j=1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
%         A(po,po+1)=9*dy2*21*cos(6*pi*(0+dy*j)); 
        A(po,po-1)=9*dy2;
%         A(po,po-(Nx-1))=4*dx2*21*cos(4*pi*(0+dx*i));
        A(po,po+(Nx-1))=4*dx2;
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-4*dx2*21*cos(4*pi*(0+dx*i))-9*dy2*21*cos(6*pi*(0+dy*j));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));

% Top right edge of the grid    
i=Nx-1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
%         A(po,po+1)=9*dy2*21*cos(6*pi*(0+dy*j)); 
        A(po,po-1)=9*dy2;
        A(po,po-(Nx-1))=4*dx2;
%         A(po,po+(Nx-1))=4*dx2*21*cos(4*pi*(0+dx*i));
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-4*dx2*21*cos(4*pi*(0+dx*i))-9*dy2*21*cos(6*pi*(0+dy*j));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));

% Top left edge of the grid        
i=1;j=Ny-1;
        po=i+(j-1)*(Nx-1);
        A(po,po)= -2*(9*dy2+4*dx2);
        A(po,po+1)=9*dy2; 
%         A(po,po-1)=9*dy2*21*cos(6*pi*(0+dy*j));
        A(po,po-(Nx-1))=4*dx2;
%         A(po,po+(Nx-1))=4*dx2*21*cos(4*pi*(0+dx*i));
        
        b(po)=dx2dy2*(-36)*21*8*(pi^2)*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)))-4*dx2*21*cos(4*pi*(0+dx*i))-9*dy2*21*cos(6*pi*(0+dy*j));
        r(po)=21*cos(2*pi*(2*(0+i*dx)+3*(0+j*dy)));
        
end