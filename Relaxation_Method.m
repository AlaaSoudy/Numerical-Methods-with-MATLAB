% Relaxation Method 

format shortg;
fprintf('\n====================  Relaxation Method  ====================\n\n');
n =input('Enter number of equations: ');
fprintf('Enter the augmented matrix [A|b] row by row:\n');
aug=zeros(n, n+1);
for i = 1:n
    aug(i,:) = input(sprintf('Row %d: ', i));
end
A= aug(:, 1:n);
b= aug(:, end);
if ~isDiagonallyDominant(A)
    fprintf('Warning: The system is NOT strictly diagonally dominant.\n');
    fprintf('Convergence is NOT guaranteed!\n\n');
else
    fprintf('Good: The system is diagonally dominant → Convergence expected.\n\n');
end
tol_input=input('Tolerance : ', 's');
if isempty(tol_input)
    tol = 1e-6;
else
    tol = str2double(tol_input);
    if isnan(tol) || tol <= 0
        tol = 1e-6;
    end
end
x_input = input('Initial guess (press Enter for zeros): ');
if isempty(x_input)
    x = zeros(n,1);
else
    x=x_input(:);
end
fprintf('  k   ');
for i = 1:n, fprintf('     x%-7d', i);
end
for i = 1:n, fprintf('     r%-7d', i); 
end
fprintf('     Error\n');
fprintf('%s\n', repmat('-', 1, 5 + 13*n*2 + 12));
k=0;
max_iter = 1000;
while true
    r_raw = b - A*x;
    r = r_raw ./ diag(A);          
    max_r = max(abs(r));
    fprintf('%3d  ', k);
    for i = 1:n,fprintf('%12.6f ', x(i)); 
    end
    for i = 1:n,fprintf('%12.6f ', r(i)); 
    end
    fprintf('%12.6f\n', max_r);
    if max_r < tol
        break;
    end
    if k >= max_iter
        fprintf('\nWarning: Maximum iterations (%d) reached!\n', max_iter);
        break;
    end
    [~, idx] = max(abs(r));
    x(idx) = x(idx) + r(idx);   
    k=k + 1;
end
printSolution(x, k,tol);

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