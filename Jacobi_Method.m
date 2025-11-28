% Jacobi Method

clc; clear;
fprintf('\n==================== Jacobi Method ====================\n\n');
n=input('Enter number of equations: ');
fprintf('Enter the augmented matrix [A|b] row by row:\n');
M=zeros(n,n+1);
for i = 1:n
    M(i,:)=input(sprintf('Row %d: ',i));
end
A=M(:,1:n); b=M(:,end);
if ~isDiagonallyDominant(A)
    fprintf('Warning: The system is NOT strictly diagonally dominant.\n');
    fprintf('Convergence is NOT guaranteed!\n\n');
else
    fprintf('Good: The system is diagonally dominant → Convergence expected.\n\n');
end
tol = input('Tolerance : ');
if isempty(tol), tol = 0.01;
end
while tol<= 0, tol = input('Enter positive tolerance: '); 
end
k=0;
x=input('Initial guess (press Enter for zeros): ');
if isempty(x)
    x = zeros(n,1); 
else 
    x = x(:); 
end
fprintf('\n');
fprintf('\n');
fprintf(' %-4s ', 'k');
for i= 1:n
   fprintf(' %-14s ', sprintf('x%d',i));
end
fprintf(' %-20s\n', 'Error');
fprintf('%s\n', repmat('-', 1, 6 + 16*n + 20));
fprintf(' %-4d ', k);
for j= 1:n
    fprintf(' %-14.8f ', x(j));
end
fprintf(' %s \n', '----');
while true
    k= k+1;
    xold=x;
    for i = 1:n
        idx= [1:i-1, i+1:n]; 
        x(i)= (b(i) - A(i,idx) * xold(idx)) / A(i,i);
       
        if abs(A(i,i)) < eps
            error('Diagonal element A(%d,%d) is zero! Jacobi method cannot proceed.', i, i);
        end
    end
    err=norm(x-xold, inf); 
   fprintf(' %-4d ', k);
for j = 1:n
    fprintf(' %-14.8f ', x(j));
end
    if err<tol
        fprintf(' %-16.8f < %.3f\n', err,tol);
        break;
    else
        fprintf(' %-16.8f \n', err);
    end
end
printSolution(x,k,tol)

% ===================== Functions =====================
function dominant=isDiagonallyDominant(A)
    n = size(A,1);
    dominant = true;
    for i = 1:n
        if abs(A(i,i)) <= sum(abs(A(i,:))) - abs(A(i,i))
            dominant = false;
            return;
        end
    end
end
function printSolution(x, k, tol)
    fprintf('\nConverged in %d iterations (tol = %.5g)\n', k, tol);
    fprintf('%s\n', repmat('=',1,50));
    for i = 1:length(x)
        fprintf('   x%d = %.10f\n', i, x(i));
    end
    fprintf('   Rounded: ');
    for i = 1:length(x)
        fprintf('x%d = %.4f   ', i, x(i));
    end
    fprintf('\b\b\b\n');
    fprintf('%s\n', repmat('=',1,50));
end