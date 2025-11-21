% Fixed Point Iteration

g= input('Enter your function g(x): ','s');
g= str2func(['@(x) ' g]); 
x0= input('Enter initial value: ');
e=input('Enter tolerence: ');
n=input('Number of itteration: ');
converged = false; 
   fprintf('%-6s %-16s %-16s %-16s\n','Iter','xᵢ₋₁','xᵢ','Error');
fprintf('%-6d %-16.10f %-16.10f %-16.10f\n',0,0,x0,abs(x0-0));
for i=1:n
    x1 =g(x0);
    err=abs(x1-x0);
    fprintf('%-6d %-16.10f %-16.10f %-16.10f\n', i,x0,x1,err);
    if err<e
       fprintf('Root found ≈ %.10f\n',x1);
       converged =true;
    break;
     end
        x0=x1;
end

if ~converged
    fprintf('\nMethod did not converge within %d iterations.\n',n);
end