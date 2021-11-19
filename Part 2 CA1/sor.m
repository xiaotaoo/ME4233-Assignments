function [residual_s, iter_s]=sor(A,b, soln,u0,w)
D=diag(diag(A)); 
L_S=tril(A,-1);     %take the lower triangular matrix of A excluding the diagonal
U_S=triu(A,1);      %take the upper triangular matrix of A excluding the diagonal
iter_s=0;           %initialize the iteration for S
residual_s=[];      %initialize the residual array
while 1
    u1=(D+w*L_S)\(w*b-(w*U_S+(w-1)*D)*u0);  %the general formula for the SOR method 
    mul = (u0-soln).*(u0-soln);%array multiplication to get the square of u(guess)-utrue(solnQR)
    res_s = sqrt(sum(mul)); %obtain the magnitude of the difference in vectors
    residual_s = [residual_s res_s];    %append the residual into the list for plotting
    if res_s<10^-7 %CHANGE TO 10^-17 IF SUBMITTING
        break
    end
    u0=u1;
    iter_s=iter_s+1;
    disp(['Iteration number:' num2str(iter_s) ' residual=' num2str(res_s)])
end