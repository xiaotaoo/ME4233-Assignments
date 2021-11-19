function [vort_2] = adv_vort(method,streammfunc,vort,Nx,Ny,dx,dy,dt,Re,t)
    %Create split for Explicit and Implict Methods
    if method == 1 %Use Explicit Euler
        %Direct Adjustment of Matrix
        RHS = assembleRHS(Nx,Ny,streammfunc,vort,Re,dx,dy,t);
        vort_2 = vort + dt*RHS;
    elseif method == 2 %Use Implicit Euler
        %Assembly Au = B for vorticity
        [A,b] = assembleAb_implicit(Nx,Ny,streammfunc,vort,Re,dx,dy,dt,t);
        Ngs=(Nx-1)*(Ny-1);
        u0=zeros(Ngs,1);
        %Inverse Matrix
        vort_2=SOR(A,b,u0);
        vort_2=reshape(vort_2,Nx-1,Ny-1);
    end                        
end









