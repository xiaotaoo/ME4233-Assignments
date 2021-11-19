function [soln_LU] = LU(A,b)

%LU Decomposition
N= size(A,1);               %obtain the integer N which corresponds to the number of rows and columns of A
L= eye(N,N);                %create an identity matrix of size NxN
A_LU = sparse(A);                   %create a separate array A_LU that will undergo the gaussian elimination

for i=1:N-1                 %perform gaussian elimination
    Lo=eye(N,N);            %create an identity Lo matrix for each iteration of the gaussian elimination
    Lo(i+1:N,i) = - A_LU(i+1:N,i)/A_LU(i,i); %Use the LU decomposition formula to fill up all elements below the diagonal for each column
    A_LU(i+1:N,:) = A_LU(i+1:N,:) + Lo(i+1:N,i)*A_LU(i,:); %Subtract the Lo from each element under the diagonal column
    L=Lo*L;              %left multiply the original identity matrix by the L matrix obtained from each iteration
    disp(num2str(i))
end

U = A_LU;

y_LU = L*b; %Obtain the intermediate vector y
soln_LU=U\y_LU;