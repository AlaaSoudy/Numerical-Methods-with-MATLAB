% Secant Method

f = input('Enter your function f(x): ','s');
f = str2func(['@(x) ' f]); 
x0= input('Enter initial value x0: ');
x1=input('Enter initial value x1: ');
e=input('Enter tolerence: ');
n=input('Number of itteration: ');
fprintf('%-6s %-16s %-16s %-16s %-16s\n','Iter','xᵢ₋2','xᵢ₋₁','xᵢ','Error');
for i= 1:n
    x2 = x1-( (x1-x0) / (f(x1)-f(x0)))*f(x1); 
     err=abs(x2-x1);
     fprintf('%-6d %-16.10f %-16.10f %-16.10f %-16.10f\n', i,x0, x1,x2, err);
  if  err< e
  fprintf('Root found ≈ %.10f\n',x2);
       break;
 end
    x0=x1;  
    x1=x2;
end