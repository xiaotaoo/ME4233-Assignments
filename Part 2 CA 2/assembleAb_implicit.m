function [A,b]=assembleAb_implicit(Nx,Ny,streamfunc,vort,Re,dx,dy,dt,t)
    % Init BC
        U_south = 0;  
        U_north = sin(2*pi*5*t);
        U_west  = 0; 
        U_east  = 0; 
    % Init constants
        nu=1/Re;
        dx2=dx*dx;
        dy2=dy*dy;
    % Init A & B matrix
        A=zeros((Nx-1)*(Ny-1),(Nx-1)*(Ny-1));
        b=zeros((Nx-1)*(Ny-1),1);
    % interior
    for j=2:Ny-2
        for i=2:Nx-2
            % streamfunction
            fac1 = -(streamfunc(i,j+1) - streamfunc(i,j-1))/2/dy;    
            fac2 =  (streamfunc(i+1,j) - streamfunc(i-1,j))/2/dx;
            % RHS at i,j
            po=i+(j-1)*(Nx-1);
            A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
            A(po,po+1)=(-fac1/2/dx)+(-nu/dx2); %i+1,j
            A(po,po-1)=(fac1/2/dx)+(-nu/dx2); %i-1,j
            A(po,po-(Nx-1))=(fac2/2/dy)+(-nu/dy2); %i,j-1
            A(po,po+(Nx-1))=(-fac2/2/dy)+(-nu/dy2); %i,j+1
            %include vort(i,j) at prev timestep
            b(po)=(1/dt)*vort(i,j);
        end
    end
    % South
    j=1;
    for i=2:Nx-2
        vortsouthbc = ( 0 - streamfunc(i,j) - U_south*dy)/(0.5*dy*dy);
        % streamfunction
        fac1 = -(streamfunc(i,j+1) - 0             )/2/dy;    
        fac2 =  (streamfunc(i+1,j) - streamfunc(i-1,j))/2/dx;
        % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po+1)=(-fac1/2/dx)+(-nu/dx2); %i+1,j
        A(po,po-1)=(fac1/2/dx)+(-nu/dx2); %i-1,j
        A(po,po+(Nx-1))=(-fac2/2/dy)+(-nu/dy2); %i,j+1
        %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) +...
            (((-fac2/2/dy)+(nu/dy2))*vortsouthbc);
    end

    % North
    j=Ny-1;
    for i=2:Nx-2
        vortnorthbc = ( 0 - streamfunc(i,j) - U_north*dy)/(0.5*dy*dy);
        % streamfunction
        fac1 = -(0              - streamfunc(i,j-1))/2/dy;    
        fac2 =  (streamfunc(i+1,j) - streamfunc(i-1,j))/2/dx;
        % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po+1)=(-fac1/2/dx)+(-nu/dx2); %i+1,j
        A(po,po-1)=(fac1/2/dx)+(-nu/dx2); %i-1,j
        A(po,po-(Nx-1))=(fac2/2/dy)+(-nu/dy2); %i,j-1
        %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) +...
            (((fac2/2/dy)+(nu/dy2))*vortnorthbc);
    end
    % West
    i=1;
    for j=2:Ny-2
        vortwestbc = ( 0 - streamfunc(i,j) - U_west*dx)/(0.5*dx*dx);
        % streamfunction
        fac1 = -(streamfunc(i,j+1) - streamfunc(i,j-1))/2/dy;    
        fac2 =  (streamfunc(i+1,j) - 0             )/2/dx;
        % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po+1)=(-fac1/2/dx)+(-nu/dx2); %i+1,j
        A(po,po-(Nx-1))=(fac2/2/dy)+(-nu/dy2); %i,j-1
        A(po,po+(Nx-1))=(-fac2/2/dy)+(-nu/dy2); %i,j+1
        %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) +...
            (((-fac1/2/dx)+(nu/dx2))*vortwestbc);
    end
    % East
    i=Nx-1;
    for j=2:Ny-2
            vorteastbc = ( 0 - streamfunc(i,j) - U_east*dx)/(0.5*dx*dx);
            % streamfunction
            fac1 = -(streamfunc(i,j+1) - streamfunc(i,j-1))/2/dy;    
            fac2 =  (0              - streamfunc(i-1,j))/2/dx;        
            % RHS at i,j
            po=i+(j-1)*(Nx-1);
            A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
            A(po,po-1)=(fac1/2/dx)+(-nu/dx2); %i-1,j
            A(po,po-(Nx-1))=(fac2/2/dy)+(-nu/dy2); %i,j-1
            A(po,po+(Nx-1))=(-fac2/2/dy)+(-nu/dy2); %i,j+1
            %include vort(i,j) at prev timestep
            b(po)=(1/dt)*vort(i,j) +...
                (((fac1/2/dx)+(nu/dx2))*vorteastbc);
    end
    % South-west
        i=1;j=1;
        vortsouthbc = ( 0 - streamfunc(i,j) - U_south*dy)/(0.5*dy*dy);
        vortwestbc = ( 0 - streamfunc(i,j) - U_west*dx)/(0.5*dx*dx);
    % streamfunction
        fac1 = -(streamfunc(i,j+1) - 0 )/2/dy;    
        fac2 =  (streamfunc(i+1,j) - 0 )/2/dx;
    % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po+1)=(-fac1/2/dx)+(-nu/dx2); %i+1,j
        A(po,po+(Nx-1))=(-fac2/2/dy)+(-nu/dy2); %i,j+1
    %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) +...
            (((-fac1/2/dx)+(nu/dx2))*vortwestbc) +...
            (((-fac2/2/dy)+(nu/dy2))*vortsouthbc);
    % South-east          
        i=Nx-1;j=1;
        vortsouthbc = ( 0 - streamfunc(i,j) - U_south*dy)/(0.5*dy*dy);
        vorteastbc = ( 0 - streamfunc(i,j) - U_east*dx)/(0.5*dx*dx);
    % streamfunction
        fac1 = -(streamfunc(i,j+1) - 0)/2/dy;    
        fac2 =  (0 - streamfunc(i-1,j))/2/dx;
    % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po-1)=(fac1/2/dx)+(-nu/dx2); %i-1,j
        A(po,po+(Nx-1))=(-fac2/2/dy)+(-nu/dy2); %i,j+1
    %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) +...
            (((-fac2/2/dy)+(nu/dy2))*vortsouthbc) +...
            (((fac1/2/dx)+(nu/dx2))*vorteastbc);
    % North-east         
        i=Nx-1;j=Ny-1;
        vortnorthbc = (0-streamfunc(i,j) - U_north*dy)/(0.5*dy2);
        vorteastbc =  (0-streamfunc(i,j) - U_east*dx)/(0.5*dx2);
    % streamfunction
        fac1 = -(0 - streamfunc(i,j-1))/2/dy;    
        fac2 =  (0 - streamfunc(i-1,j))/2/dx;
    % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po-1)=(fac1/2/dx)+(-nu/dx2); %i-1,j
        A(po,po-(Nx-1))=(fac2/2/dy)+(-nu/dy2); %i,j-1
    %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) + ...
            (((fac2/2/dy)+(nu/dy2))*vortnorthbc) + ...
            (((fac1/2/dx)+(nu/dx2))*vorteastbc);
    % North-west          
        i=1;j=Ny-1;
        vortnorthbc = ( 0 - streamfunc(i,j) - U_north*dy)/(0.5*dy*dy);
        vortwestbc  = ( 0 - streamfunc(i,j) + U_west *dx)/(0.5*dx*dx);
    % streamfunction
        fac1 = -(0 - streamfunc(i,j-1))/2/dy;    
        fac2 =  (streamfunc(i+1,j) - 0)/2/dx;
    % RHS at i,j
        po=i+(j-1)*(Nx-1);
        A(po,po)= (1/dt) + (2*nu/dx2) + (2*nu/dy2);
        A(po,po+1)=(-fac1/2/dx)+(-nu/dx2); %i+1,j
        A(po,po-(Nx-1))=(fac2/2/dy)+(-nu/dy2); %i,j-1
    %include vort(i,j) at prev timestep
        b(po)=(1/dt)*vort(i,j) + ...
            (((fac2/2/dy)+(nu/dy2))*vortnorthbc) + ...
            (((-fac1/2/dx)+(nu/dx2))*vortwestbc);
end