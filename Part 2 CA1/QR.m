function [soln_QR] = QR(A,b)
N= size(A,1);               %obtain the integer N which corresponds to the number of rows and columns of A
Q= eye(N,N);
Q(:,1)=A(:,1);  %Equate the first column of vector Q = first column of vector A
for j=2:N
    P=zeros(N,1);   %define an empty column vector every loop
    for k=1:j-1     %loop through every column before the J column (j-1)
        P=P+(Q(:,k)'*A(:,j))/(Q(:,k)'*Q(:,k))*Q(:,k);   %we are projecting the first to j-1 column of Q on the j column of A then summing up the projections
    end
    Q(:,j)=A(:,j)-P; %subtract the sum of the projections from the j column of A, obtaining v1
    Q=sparse(Q);
    disp(['j:' num2str(j)])
end
for col=1:N
    Q(:,col)=Q(:, col)/sqrt(Q(:,col)'*Q(:,col)); %loop through the columns in Q to normalize the new vectors, obtaining the e matrix
        disp(['col:' num2str(col)])
end

R = zeros(N,N);
for c=1:N
    for d=1:c %loop from 1 to the current column index, c
        R(d,c)=A(:,c)'*Q(:,d); %for each column in the R matrix, add the <u,e> into the corresponding elements above the diagonal
    end
    disp(['c:' num2str(c)])
end         %We can equate A=QR after obtaining the Q and R matrix
y_QR=Q' *b; %Transpose the Q matrix since the inverse of Q is just the transpose and multiply it by B, Ru=QT b
soln_QR=R\y_QR;  %u= R(-1) * QT * B