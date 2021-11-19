function [RHS] = assembleRHS(Nx,Ny,streamfunc,vort,Re,dx,dy,t)
    % Init BC
        U_south = 0;  
        U_north = sin(2*pi*5*t);
        U_west  = 0; 
        U_east  = 0; 
    % Init Constants and RHS Matrix
        nu=1/Re;
        RHS=zeros(Nx-1,Ny-1);
    % interior
    for j=2:Ny-2
        for i=2:Nx-2
            fac1 = -(streamfunc(i,j+1) - streamfunc(i,j-1))/2/dy;
            fac2 =  (streamfunc(i+1,j) - streamfunc(i-1,j))/2/dx;
            RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
                        fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
                        nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/dx/dx + ... 
                             ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/dy/dy );                     
        end
    end
    % south
    j=1;
    for i=2:Nx-2
        vortsouthbc = ( 0 - streamfunc(i,j) + U_south*dy)/(0.5*dy*dy);    
        fac1 = -(streamfunc(i,j+1) - 0             )/2/dy;
        fac2 =  (streamfunc(i+1,j) - streamfunc(i-1,j))/2/dx;
        RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
                    fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
                    nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/dx/dx + ... 
                         ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/dy/dy ); 
    end
    % North
    j=Ny-1;
    for i=2:Nx-2
        vortnorthbc = ( 0 - streamfunc(i,j) - U_north*dy)/(0.5*dy*dy);
        fac1 = -(0              - streamfunc(i,j-1))/2/dy;
        fac2 =  (streamfunc(i+1,j) - streamfunc(i-1,j))/2/dx;
        RHS(i,j) =  fac1*( vort(i+1,j)-vort(i-1,j) )/2/dx + ...
                    fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
                    nu*( ( vort(i+1,j)-2*vort(i,j)+vort(i-1,j) )/dx/dx + ... 
                         ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/dy/dy ); 
    end
    % West
    i=1;
    for j=2:Ny-2
        vortwestbc = ( 0 - streamfunc(i,j) - U_west*dx)/(0.5*dx*dx);
        fac1 = -(streamfunc(i,j+1) - streamfunc(i,j-1))/2/dy;
        fac2 =  (streamfunc(i+1,j) - 0             )/2/dx;
        RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
                    fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
                    nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/dx/dx + ... 
                         ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/dy/dy ); 
    end
    % East
    i=Nx-1;
    for j=2:Ny-2
        vorteastbc = ( 0 - streamfunc(i,j) + U_east*dx)/(0.5*dx*dx);
        fac1 = -(streamfunc(i,j+1) - streamfunc(i,j-1))/2/dy;
        fac2 =  (0              - streamfunc(i-1,j))/2/dx;
        RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
                    fac2*( vort(i,j+1)-vort(i,j-1) )/2/dy + ...
                    nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/dx/dx + ... 
                         ( vort(i,j+1)-2*vort(i,j)+vort(i,j-1) )/dy/dy ); 
    end
    % South-west
        i=1;j=1;
        vortsouthbc = ( 0 - streamfunc(i,j) + U_south*dy)/(0.5*dy*dy);
        vortwestbc  = ( 0 - streamfunc(i,j) + U_west *dx)/(0.5*dx*dx);
        fac1 = -(streamfunc(i,j+1) - 0 )/2/dy;
        fac2 =  (streamfunc(i+1,j) - 0 )/2/dx;
        RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
                    fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
                    nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/dx/dx + ... 
                         ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/dy/dy ); 

    % South-east          
        i=Nx-1;j=1;
        vortsouthbc = ( 0 - streamfunc(i,j) + U_south*dy)/(0.5*dy*dy);
        vorteastbc  = ( 0 - streamfunc(i,j) - U_east *dx)/(0.5*dx*dx);
        fac1 = -(streamfunc(i,j+1) - 0)/2/dy;
        fac2 =  (0 - streamfunc(i-1,j))/2/dx;
        RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
                    fac2*( vort(i,j+1)-vortsouthbc )/2/dy + ...
                    nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/dx/dx + ... 
                         ( vort(i,j+1)-2*vort(i,j)+vortsouthbc )/dy/dy );              
    % North-east         
        i=Nx-1;j=Ny-1;
        vortnorthbc = ( 0 - streamfunc(i,j) - U_north*dy)/(0.5*dy*dy);
        vorteastbc  = ( 0 - streamfunc(i,j) - U_east *dx)/(0.5*dx*dx); 
        fac1 = -(0 - streamfunc(i,j-1))/2/dy;
        fac2 =  (0 - streamfunc(i-1,j))/2/dx;
        RHS(i,j) =  fac1*( vorteastbc -vort(i-1,j) )/2/dx + ...
                    fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
                    nu*( ( vorteastbc -2*vort(i,j)+vort(i-1,j) )/dx/dx + ... 
                         ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/dy/dy );
    % North-west          
        i=1;j=Ny-1;
        vortnorthbc = ( 0 - streamfunc(i,j) - U_north*dy)/(0.5*dy*dy);
        vortwestbc  = ( 0 - streamfunc(i,j) + U_west *dx)/(0.5*dx*dx);
        fac1 = -(0 - streamfunc(i,j-1))/2/dy;
        fac2 =  (streamfunc(i+1,j) - 0)/2/dx;
        RHS(i,j) =  fac1*( vort(i+1,j)-vortwestbc  )/2/dx + ...
                    fac2*( vortnorthbc-vort(i,j-1) )/2/dy + ...
                    nu*( ( vort(i+1,j)-2*vort(i,j)+vortwestbc  )/dx/dx + ... 
                         ( vortnorthbc-2*vort(i,j)+vort(i,j-1) )/dy/dy );                                     
end

