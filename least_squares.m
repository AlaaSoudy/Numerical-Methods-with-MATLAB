%CURVE FITTING - LEAST SQUARES
clc; clear; close all;
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
        E = sprintf('y=%.6f*x^{%.6f}',A,B);
    else
        E = sprintf('y=%.6f*x^{-%.6f}',A,abs(B));
    end
    YF = A*X.^B;
    R2 = calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_exp(X,Y)
    if any(Y<=0)
        error('Exponential fit needs positive y!');
    end
    [lnA,B] = lin_fit(X,log(Y));
    A = exp(lnA);
    C = [A,B];
    if B>=0
        E = sprintf('y=%.6f*exp(%.6f*x)',A,B);
    else
        E = sprintf('y=%.6f*exp(-%.6f*x)',A,abs(B));
    end
    YF = A*exp(B*X);
    R2 = calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_rat(X,Y)
    if any(Y==0)
        error('Rational fit needs non-zero y!');
    end
    [A,B] = lin_fit(X,X./Y);
    C = [A,B];
    if B>=0
        E = sprintf('y=x/(%.6f+%.6f*x)',A,B);
    else
        E = sprintf('y=x/(%.6f-%.6f*x)',A,abs(B));
    end
    YF = X./(A+B.*X);
    R2 = calc_R2(Y,YF);
end
function [C,E,R2,YF] = fit_cust(X,Y,eq)
    prm = {'a','b','c','d','e','f'};
    idx = [];
    for i=1:length(prm)
        if contains(eq,prm{i})
            idx = [idx,i];
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
            rep = sprintf('(%.6f)',val);
        else
            rep = sprintf('%.6f',val);
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
    fs = strrep(fs,'^','.^');
    fs = strrep(fs,'*','.*');
    fs = strrep(fs,'/','./');
    fh = str2func(['@(p,x)' fs]);
end
function R2 = calc_R2(Ya,Yp)
    SSr = sum((Ya-Yp).^2);
    SSt = sum((Ya-mean(Ya)).^2);
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
  
end