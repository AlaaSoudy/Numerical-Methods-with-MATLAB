%NEWTON DIVIDED DIFFERENCES COMPARISON
clc; clear; close all;
disp('===========================================');
disp('NEWTON DIVIDED DIFFERENCES COMPARISON');
disp('===========================================');
n = input('Enter number of points: ');
x = input('Enter x values (e.g., [0 1 2 5]): ');
y = input('Enter y values (e.g., [2 3 12 147]): ');
if length(x) ~= n || length(y) ~= n
    error('Number of points does not match input vectors!');
end
[x, y] = validate_data(x, y);
display_points(x, y, 'SORTED DATA POINTS');
[coeffs_fwd, poly_fwd, table_fwd] = newton_forward_dd(x, y);
[coeffs_bwd, poly_bwd, table_bwd] = newton_backward_dd(x, y);
fprintf('\n===========================================\n');
fprintf('COMPARISON OF NEWTON METHODS\n');
fprintf('===========================================\n');
fprintf('\n1. NEWTON FORWARD:\n');
print_dd_table(table_fwd, x, y, 'forward');
fprintf('\nCoefficients:\n');
for i = 1:length(coeffs_fwd)
    fprintf('f[x0');
    for j = 1:i-1
        fprintf(',x%d', j);
    end
    fprintf('] = %.6f\n', coeffs_fwd(i));
end
fprintf('\nPolynomial: %s\n', poly_fwd);

fprintf('\n2. NEWTON BACKWARD:\n');
print_dd_table(table_bwd, x, y, 'backward');
fprintf('\nCoefficients:\n');
for i = 1:length(coeffs_bwd)
    fprintf('f[x%d', n-i);
    for j = n-i+1:n-1
        fprintf(',x%d', j);
    end
    fprintf('] = %.6f\n', coeffs_bwd(i));
end
fprintf('\nPolynomial: %s\n', poly_bwd);

fprintf('\n===========================================\n');
fprintf('INTERPOLATION\n');
fprintf('===========================================\n');

while true
    x_input = input('\nEnter x value (or Enter to exit): ', 's');
    if isempty(x_input)
        fprintf('\nProgram completed.\n');
        break;
    end
    
    x_val = str2double(x_input);
    if isnan(x_val)
        fprintf('Invalid input! Enter numeric value.\n');
        continue;
    end
    
    y_fwd = evaluate_forward_dd(x_val, x, coeffs_fwd);
    y_bwd = evaluate_backward_dd(x_val, x, coeffs_bwd);
    
    fprintf('\nFor x = %.4f:\n', x_val);
    fprintf(' Forward:  y = %.6f\n', y_fwd);
    fprintf(' Backward: y = %.6f\n', y_bwd);
    
    diff_val = abs(y_fwd - y_bwd);
    fprintf(' Difference: %.2e\n', diff_val);
    
    if diff_val > 1e-10
        fprintf(' Warning: Large difference detected!\n');
    end
end
function [coeffs, poly_str, dd_table] = newton_forward_dd(x, y)
    n = length(x);
    dd_table = zeros(n, n);
    dd_table(:, 1) = y(:);
        for j = 2:n
        for i = 1:n-j+1
            dd_table(i, j) = (dd_table(i+1, j-1) - dd_table(i, j-1)) / (x(i+j-1) - x(i));
        end
    end
        coeffs = dd_table(1, :);
    poly_str = build_polynomial(x, coeffs, 'forward');
end
function [coeffs, poly_str, dd_table] = newton_backward_dd(x, y)
    n = length(x);
    dd_table = zeros(n, n);
    dd_table(:, 1) = y(:);
    
    for j = 2:n
        for i = 1:n-j+1
            dd_table(i, j) = (dd_table(i+1, j-1) - dd_table(i, j-1)) / (x(i+j-1) - x(i));
        end
    end
        coeffs = zeros(1, n);
    for i = 1:n
        coeffs(i) = dd_table(n-i+1, i);
    end
    
    poly_str = build_polynomial(x, coeffs, 'backward');
end
function poly_str = build_polynomial(x, coeffs, method)
    n = length(coeffs);
    poly_str = 'P(x) = ';
    first_term = true;
        for i = 1:n
        coeff = coeffs(i);
        if abs(coeff) > 1e-12
                 coeff_abs = abs(coeff);
                  if abs(coeff_abs - 1) < 1e-10
                coeff_str = '';
                if first_term && coeff > 0
                    coeff_str = '1';
                end
            else
                coeff_str = sprintf('%.6f', coeff_abs);
            end
                    if ~first_term
                if coeff > 0
                    poly_str = [poly_str ' + '];
                else
                    poly_str = [poly_str ' - '];
                end
            elseif coeff < 0
                poly_str = [poly_str '-'];
            end
                     if ~isempty(coeff_str)
                poly_str = [poly_str coeff_str];
            end
                        for j = 1:i-1
                if method(1) == 'f'
                    xi = x(j);
                else
                    xi = x(n-j+1);
                end
                
                if isempty(coeff_str) && j == 1
                    poly_str = [poly_str '(x-' sprintf('%.4f', xi) ')'];
                else
                    poly_str = [poly_str '*(x-' sprintf('%.4f', xi) ')'];
                end
            end
                        first_term = false;
        end
    end
    
    if first_term
        poly_str = [poly_str '0'];
    end
end
function print_dd_table(table, x, y, method)
    n = length(x);
        if method(1) == 'f'
        symbol = 'Δ';
        fprintf('\n   i       xi        f(xi)');
        for j = 1:n-1
            fprintf('      %s^%df', symbol, j);
        end
        fprintf('\n');
        fprintf('   %s\n', repmat('-', 1, 25 + (n-1)*12));
        
        for i = 1:n
            fprintf('   %2d   %8.4f   %10.6f', i, x(i), y(i));
            for j = 2:n-i+1
                fprintf('   %10.6f', table(i, j));
            end
            fprintf('\n');
        end
    else
        symbol = '∇';
        fprintf('\n   i       xi        f(xi)');
        for j = 1:n-1
            fprintf('      %s^%df', symbol, j);
        end
        fprintf('\n');
        fprintf('   %s\n', repmat('-', 1, 25 + (n-1)*12));
        
        for i = 1:n
            fprintf('   %2d   %8.4f   %10.6f', i, x(i), y(i));
   for j = 1:n-1
                if j < i
                    fprintf('   %10.6f', table(i-j, j+1));
                else
                    fprintf('   %10s', ' ');
                end
            end
            fprintf('\n');
        end
    end
end
function y_val = evaluate_forward_dd(x_val, x, coeffs)
    n = length(coeffs);
    y_val = coeffs(1);
    product = 1;
    
    for i = 2:n
        product = product * (x_val - x(i-1));
        y_val = y_val + coeffs(i) * product;
    end
end
function y_val = evaluate_backward_dd(x_val, x, coeffs)
    n = length(coeffs);
    y_val = coeffs(1);
    product = 1;
    
    for i = 2:n
        product = product * (x_val - x(end-i+2));
        y_val = y_val + coeffs(i) * product;
    end
end
function [x_valid, y_valid] = validate_data(x, y)
    if length(x) ~= length(y)
        error('x and y must have the same length');
    end
    if length(unique(x)) ~= length(x)
        error('x values must be unique');
    end
    if ~issorted(x)
        [x_valid, sort_idx] = sort(x);
        y_valid = y(sort_idx);
        fprintf('Data sorted in ascending order.\n');
    else
        x_valid = x;
        y_valid = y;
    end
end
function display_points(x, y, title_str)
    n = length(x);
    fprintf('\n%s\n', title_str);
    fprintf('%s\n', repmat('=', length(title_str), 1));
    fprintf('Number of points: %d\n\n', n);
    fprintf('   i     x(i)        y(i)   \n');
    fprintf('  --------------------------\n');
    for i = 1:n
        fprintf('  %2d   %8.4f    %8.4f\n', i, x(i), y(i));
    end

end