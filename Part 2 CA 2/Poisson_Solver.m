function [streamfunc] = Poisson_Solver(vort,A,Nx,Ny)
    b = -vort(:);
    Ngs=(Nx-1)*(Ny-1);
    u0=zeros(Ngs,1); %Init U0
    A=sparse(A);
    streamfunc = SOR(A,b,u0);
    streamfunc = reshape(streamfunc,Nx-1,Ny-1);
end

