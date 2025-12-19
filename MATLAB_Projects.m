%LAGRANGE INTERPOLATION

clc; clear; close all;
disp('===========================================');
disp('LAGRANGE INTERPOLATION');
disp('===========================================');
n=input('Enter number of points: ');
disp('Enter x values as vector : ');
x=input('x=');
if length(x)~=n
    error('Number of points does not match x vector length!');
end
disp('Enter y values as vector : ');
y=input('y=');
if length(y)~=n
    error('Number of points does not match y vector length!');
end
[x, y] = validate_data_L(x, y);
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
    if abs(denominator)<1e-12
        basis_polys{i} ='0';
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
function [x_valid,y_valid]=validate_data_L(x,y)
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

%NEWTON DIVIDED DIFFERENCES COMPARISON

disp('===========================================');
disp('NEWTON DIVIDED DIFFERENCES COMPARISON');
disp('===========================================');
n = input('Enter number of points: ');
x = input('Enter x values : ');
y = input('Enter y values : ');
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
    fprintf(' Forward:  y= %.6f\n', y_fwd);
    fprintf(' Backward: y=%.6f\n', y_bwd);
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
        for j=2:n
        for i=1:n-j+1
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
        for j= 2:n
        for i= 1:n-j+1
            dd_table(i, j) = (dd_table(i+1, j-1) - dd_table(i, j-1)) / (x(i+j-1) - x(i));
        end
    end
        coeffs = zeros(1, n);
    for i = 1:n
        coeffs(i) = dd_table(n-i+1, i);
    end
        poly_str=build_polynomial(x, coeffs, 'backward');
end
function poly_str =build_polynomial(x, coeffs, method)
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
                    coeff_str='1';
                end
            else
                coeff_str =sprintf('%.6f', coeff_abs);
            end
           if ~first_term
                if coeff>0
                    poly_str=[poly_str ' + '];
                else
                    poly_str=[poly_str ' - '];
                end
            elseif coeff < 0
                poly_str=[poly_str '-'];
            end
                     if ~isempty(coeff_str)
                poly_str =[poly_str coeff_str];
            end
                        for j = 1:i-1
                if method(1) == 'f'
                    xi = x(j);
                else
                    xi = x(n-j+1);
                end
           if isempty(coeff_str) && j == 1
                    poly_str= [poly_str '(x-' sprintf('%.4f', xi) ')'];
                else
                    poly_str= [poly_str '*(x-' sprintf('%.4f', xi) ')'];
           end
            end
                first_term =false;
        end
    end
        if first_term
        poly_str = [poly_str '0'];
    end
end
function print_dd_table(table, x, y, method)
    n = length(x);
        if method(1) =='f'
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

%CURVE FITTING - LEAST SQUARES

disp('===========================================');
disp('CURVE FITTING - LEAST SQUARES');
disp('===========================================');
np = input('Enter number of points: ');
X = input('Enter x values: ');
Y = input('Enter y values: ');
if length(X)~=np || length(Y)~=np
    error('Input vectors length mismatch!');
end
display_points(X,Y,'DATA POINTS');
fprintf('\nFITTING METHODS:\n');
disp('1. Linear: y = a + b*x');
disp('2. Quadratic: y = a + b*x + c*x^2');
disp('3. Power: y = a*x^b');
disp('4. Exponential: y = a*exp(b*x)');
disp('5. Rational: y = x/(a + b*x)');
disp('6. Custom Equation');
valid = false;
while ~valid
    m = input('\nSelect method (1-6): ','s');
    met = str2double(m);
    if ~isnan(met) && isscalar(met) && met>=1 && met<=6 && mod(met,1)==0
        valid = true;
    else
        fprintf('Invalid! Enter 1-6.\n');
    end
end
fprintf('\nRESULTS:\n');
switch met
    case 1
        [C,E,R2,YF] = fit_lin(X,Y);
    case 2
        [C,E,R2,YF] = fit_quad(X,Y);
    case 3
        [C,E,R2,YF] = fit_pow(X,Y);
    case 4
        [C,E,R2,YF] = fit_exp(X,Y);
    case 5
        [C,E,R2,YF] = fit_rat(X,Y);
    case 6
        eq = input('Enter custom equation: ','s');
        [C,E,R2,YF] = fit_cust(X,Y,eq);
end
fprintf('\nEquation: %s\n',E);
fprintf('\nCoefficients:\n');
for i=1:length(C)
    fprintf('  P%d: %.6f\n',i,C(i));
end
function [A,B] = lin_fit(X,Y)
    n = length(X);
    Sx = sum(X); Sy = sum(Y); Sxx = sum(X.^2); Sxy = sum(X.*Y);
    B = (n*Sxy - Sx*Sy)/(n*Sxx - Sx^2);
    A = (Sy - B*Sx)/n;
end
function [C,E,R2,YF] = fit_lin(X,Y)
    [A,B] = lin_fit(X,Y);
    C = [A,B];
    if B>=0
        E = sprintf('y=%.6f+%.6f*x',A,B);
    else
        E = sprintf('y=%.6f-%.6f*x',A,abs(B));
    end
    YF = A + B*X;
    R2 = calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_quad(X,Y)
    n = length(X);
    M = [n,sum(X),sum(X.^2);sum(X),sum(X.^2),sum(X.^3);sum(X.^2),sum(X.^3),sum(X.^4)];
    V = [sum(Y);sum(X.*Y);sum(X.^2.*Y)];
    P = M\V;
    C = P;
        E = 'y=';
    if abs(P(1))>1e-10
        E = sprintf('%s%.6f',E,P(1));
    end
    if abs(P(2))>1e-10
        if P(2)>=0
            E = sprintf('%s+%.6f*x',E,P(2));
        else
            E = sprintf('%s-%.6f*x',E,abs(P(2)));
        end
    end
    if abs(P(3))>1e-10
        if P(3)>=0
            E = sprintf('%s+%.6f*x^2',E,P(3));
        else
            E = sprintf('%s-%.6f*x^2',E,abs(P(3)));
        end
    end
    if strcmp(E,'y=')
        E = 'y=0';
    end
    
    YF = polyval(flip(P),X);
    R2 = calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_pow(X,Y)
    if any(X<=0)||any(Y<=0)
        error('Power fit needs positive values!');
    end
    [lnA,B] = lin_fit(log(X),log(Y));
    A = exp(lnA);
    C = [A,B];
    if B>=0
        E=sprintf('y=%.6f*x^{%.6f}',A,B);
    else
        E=sprintf('y=%.6f*x^{-%.6f}',A,abs(B));
    end
    YF=A*X.^B;
    R2=calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_exp(X,Y)
    if any(Y<=0)
        error('Exponential fit needs positive y!');
    end
    [lnA,B] = lin_fit(X,log(Y));
    A=exp(lnA);
    C=[A,B];
    if B>=0
        E=sprintf('y=%.6f*exp(%.6f*x)',A,B);
    else
        E=sprintf('y=%.6f*exp(-%.6f*x)',A,abs(B));
    end
    YF=A*exp(B*X);
    R2=calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_rat(X,Y)
    if any(Y==0)
          error('Rational fit needs non-zero y!');
    end
    [A,B] = lin_fit(X,X./Y);
    C = [A,B];
    if B>=0
        E=sprintf('y=x/(%.6f+%.6f*x)',A,B);
    else
        E=sprintf('y=x/(%.6f-%.6f*x)',A,abs(B));
    end
    YF =X./(A+B.*X);
    R2=calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_cust(X,Y,eq)
    prm = {'a','b','c','d','e','f'};
    idx = [];
    for i=1:length(prm)
        if contains(eq,prm{i})
              idx=[idx,i];
        end
    end
        if isempty(idx)
        error('No parameters found!');
        end
        np = length(idx);
    p0 = ones(1,np);
    fh = make_func(eq,idx);
    opt = optimoptions('lsqcurvefit','Display','off');
    C = lsqcurvefit(fh,p0,X,Y,[],[],opt);
        E = eq;
    for i=1:np
        val = C(i);
        if val<0
            rep =sprintf('(%.6f)',val);
        else
            rep= sprintf('%.6f',val);
        end
        E = strrep(E,prm{idx(i)},rep);
    end
       YF = fh(C,X);
    R2 = calc_R2(Y,YF);
end
function fh = make_func(eq,idx)
    fs = eq;
    prm = {'a','b','c','d','e','f'};
    for i=1:length(idx)
        fs = strrep(fs,prm{idx(i)},sprintf('p(%d)',i));
    end
    fs= strrep(fs,'^','.^');
    fs=strrep(fs,'*','.*');
    fs =strrep(fs,'/','./');
    fh=str2func(['@(p,x)' fs]);
end
function R2 = calc_R2(Ya,Yp)
    SSr =sum((Ya-Yp).^2);
    SSt =sum((Ya-mean(Ya)).^2);
    if SSt==0
        R2=1;
    else
        R2=1-SSr/SSt;
    end
end
 fprintf('\nProgram completed.\n');
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
   fprintf('\n');
end