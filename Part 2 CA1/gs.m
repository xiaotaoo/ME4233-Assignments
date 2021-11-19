function [residual_g, iter_g]=gs(A,b, soln,u0)
L=tril(A);          %take the lower triangular matrix of A
U_GS = A-L;         %any remaining elements are given to the U matrix
iter_g=0;           %count the iteration for the GS method
residual_g=[];      %initialize empty array to obtain residual for GS method
while 1
    u1=L\(b-U_GS*u0);   %use the GS formula to calculate the u1 using the previous iteration
    mul = (u0-soln).*(u0-soln);%array multiplication to get the square of u(guess)-utrue(solnQR)
    res_g = sqrt(sum(mul)); %obtain the magnitude of the difference in vectors
    residual_g = [residual_g res_g];    %append the residual into the list for plotting
    if res_g<10^-7   
        break
    end
    u0=u1;      %replace the u0 with the new u1 for the next iteration
    iter_g = iter_g+1;  %update the iteration counter
    disp(['Iteration number:' num2str(iter_g) ' residual=' num2str(res_g)])
end