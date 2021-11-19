function [u0]=SOR(A,b,u0)
    w=1.5;         % the omega parameter in SOR
    resobj=10^-4;      % Residual Limit
    L = tril(A,-1);    % Get Lower Tri A
    D=diag(diag(A));   % Get the Diag A
    U=A-L-D;           % Get the Upper Tri A
    k=0;               % Iteration No.
    while 1
        u1=(D+w*L)\(w*b-(w*U+(w-1)*D)*u0);  % SOR method
        res=norm(u0-u1);            % calculate the residual through magnitude of u0 - u1
        if res<resobj              % exiting condition
            break
        end
        u0=u1;                      % updating u0 with u1
        k=k+1;                      % iteration increment
    end
end
