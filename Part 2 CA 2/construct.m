function [u,v] = construct(streamfunc,Nx,Ny,dx,dy,t)
    u =  (streamfunc(2:end-1,3:end) - streamfunc(2:end-1,1:end-2))/2/dy;
    v = -(streamfunc(3:end,2:end-1) - streamfunc(1:end-2,2:end-1))/2/dx;

    u=[0             zeros(1,Ny-1)  sin(2*pi*5*t)
       zeros(Nx-1,1) u             (sin(2*pi*5*t)*ones(Nx-1,1))
       0             zeros(1,Ny-1)  sin(2*pi*5*t) ]; 

    v=[0             zeros(1,Ny-1)  0
       zeros(Nx-1,1) v              zeros(Nx-1,1)
       0             zeros(1,Ny-1)  0 ]; 
end

