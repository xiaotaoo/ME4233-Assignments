function [residual_j, iter_j]=jacobi(A,b, soln,u0)
D=diag(diag(A));    %We have to diag twice here in order to form a diagonal matrix using the diagonal elements on A
R=A-D;              %Delete the diagonal elements from matrix A to get R
iter_j=0;             %iteration counter set to 0
residual_j=[];      %initialize an empty array to store the residual values (It will be the same number as iterations)
while 1
    u1=D\(b-R*u0);  %calculate u1 using the initial guess u0: u(k+1) = D-1 (b-Ru(k))
    mul = (u0-soln).*(u0-soln); %array multiplication to get the square of u(guess)-utrue(solnQR)
    res_j = sqrt(sum(mul)); %calculate the magnitude of the difference between the initialized vector and the calculated
    residual_j=[residual_j res_j];    %append each residual into the list for plotting
    if res_j<10^-7
        break
    end
    u0= u1;         %replace the u0 with u1 for the next iteration
    iter_j= iter_j+1;   %count the number of iterations
    disp(['Iteration number:' num2str(iter_j) ' residual=' num2str(res_j)])
end