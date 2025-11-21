% Newton-Raphson Method

syms x
f = input('Enter your function f(x): ','s');
f_sym = str2sym(f);  
df_sym = diff(f_sym, x);           
f = matlabFunction(f_sym);          
df = matlabFunction(df_sym);   
x0= input('Enter initial value: ');
e=input('Enter tolerance: ');
n=input('Number of iterations: ');
if df(x0)~=0
       fprintf('%-6s %-16s %-16s %-16s\n','Iter','xᵢ₋₁','xᵢ','Error');
    for i=1:n
        x1=x0-f(x0)/df(x0);
             err=abs(x1-x0);
       fprintf('%-6d %-16.10f %-16.10f %-16.10f\n',i,x0,x1,err);
   
      if err<e
         fprintf('Root found ≈ %.10f\n',x1);
          break;
      end
    x0=x1;
    end
else
     disp('No root.')

end
