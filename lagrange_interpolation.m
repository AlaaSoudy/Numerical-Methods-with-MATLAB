%LAGRANGE INTERPOLATION
clc; clear; close all;
disp('===========================================');
disp('LAGRANGE INTERPOLATION');
disp('===========================================');
n=input('Enter number of points: ');
disp('Enter x values as vector (e.g., [0 1 2 5]): ');
x=input('x=');
if length(x)~=n
    error('Number of points does not match x vector length!');
end
disp('Enter y values as vector (e.g., [2 3 12 147]): ');
y=input('y=');
if length(y)~=n
    error('Number of points does not match y vector length!');
end
[x, y] = validate_data(x, y);
display_points(x, y, 'INPUT DATA POINTS');
fprintf('\n===========================================\n');
fprintf('LAGRANGE INTERPOLATION RESULTS\n');
fprintf('===========================================\n');
[poly_str,basis_polys]=lagrange_interpolation(x,y);
fprintf('\nLagrange Basis Polynomials:\n');
fprintf('---------------------------\n');
for i=1:length(basis_polys)
    fprintf('L%d(x) = %s\n',i-1,basis_polys{i});
end
fprintf('\nLagrange Interpolating Polynomial:\n');
fprintf('----------------------------------\n');
fprintf('%s\n',poly_str);
fprintf('\n===========================================\n');
fprintf('INTERPOLATION TEST\n');
fprintf('===========================================\n');
while true
    fprintf('\nEnter x value for interpolation (or press Enter to exit): ');
    x_input=input('','s');
        if isempty(x_input)
        fprintf('\nProgram finished.\n');
        break;
    end
    [x_val,status]=str2num(x_input);
    if ~status
        fprintf('Invalid input! Please enter a numeric value.\n');
        continue;
    end
        y_val=evaluate_lagrange(x_val,x,y);
    fprintf('\nFor x = %.4f:\n',x_val);
    fprintf(' Lagrange Interpolation: y = %.6f\n',y_val);
end
function [poly_str,basis_polys]=lagrange_interpolation(x,y)
n=length(x);
basis_polys=cell(1,n);
for i=1:n
    numerator_parts={};
    denominator=1;
    for j=1:n
        if j~=i
            numerator_parts{end+1}=sprintf('(x-%.4f)',x(j));
            denominator=denominator*(x(i)-x(j));
        end
    end
    basis_str='';
    for k=1:length(numerator_parts)
        basis_str=[basis_str numerator_parts{k}];
        if k<length(numerator_parts)
            basis_str=[basis_str '*'];
        end
    end
    if abs(denominator) < 1e-12
        basis_polys{i} = '0';
    else
        basis_polys{i}=sprintf('(%s)/%.6f',basis_str,denominator);
    end
end
poly_str='P(x)=';
first_term=true;
for i=1:n
    if ~first_term
        if y(i) >= 0
            poly_str = [poly_str ' + '];
        else
            poly_str = [poly_str ' - '];
        end
    elseif y(i) < 0
        poly_str = [poly_str '-'];
    end
    if abs(y(i)) > 1e-10
        coeff_abs = abs(y(i));
        if abs(coeff_abs - 1) < 1e-10
            if ~first_term || y(i) > 0
                coeff_str = '';
            else
                coeff_str = '1';
            end
        else
            coeff_str = sprintf('%.6f', coeff_abs);
        end
        
        if ~isempty(coeff_str)
            poly_str = [poly_str coeff_str '*'];
        end
        poly_str = [poly_str 'L' num2str(i-1) '(x)'];
    else
        poly_str = [poly_str '0'];
    end

    first_term = false;
end
end  
function y_val=evaluate_lagrange(x_val,x,y)
n=length(x);
result=0;
for i=1:n
    term=y(i);
    for j=1:n
        if j~=i
            term=term*(x_val-x(j))/(x(i)-x(j));
        end
    end
    result=result+term;
end
y_val=result;
end
function R2=calculate_R2(y_actual,y_predicted)
SS_res=sum((y_actual-y_predicted).^2);
SS_tot=sum((y_actual-mean(y_actual)).^2);
if SS_tot==0
    R2=1;
else
    R2=1-SS_res/SS_tot;
end
end
function [a,b]=linear_fit_coeffs(x,y)
n=length(x);
Sx=sum(x);
Sy=sum(y);
Sxx=sum(x.^2);
Sxy=sum(x.*y);
b=(n*Sxy-Sx*Sy)/(n*Sxx-Sx^2);
a=(Sy-b*Sx)/n;
end
function result=factorial_custom(n)
if n<0
    error('Factorial is not defined for negative numbers');
elseif n==0||n==1
    result=1;
else
    result=1;
    for k=2:n
        result=result*k;
    end
end
end
function [x_valid,y_valid]=validate_data(x,y)
if length(x)~=length(y)
    error('x and y must have the same length');
end
if length(unique(x))~=length(x)
    error('x values must be unique');
end
if ~issorted(x)
    [x_valid,sort_idx]=sort(x);
    y_valid=y(sort_idx);
    fprintf('Note: Data has been sorted in ascending order of x.\n');
else
    x_valid=x;
    y_valid=y;
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