% False Position Method

f = input('Enter your function f(x): ','s');
f = str2func(['@(x) ' f]); 
a= input('Enter a: ');
b= input('Enter b: ');
e=input('Enter tolerance: ');
n=input('Number of iterations: ');
if a > b
    temp=a;
    a =b;
    b=temp;
end
 if f(a)*f(b)<0  
fprintf('%-6s %-14s %-14s %-14s %-16s %-16s\n','Iter','a','b','c','f(c)','Error');
   c = (a*f(b) - b*f(a)) / (f(b) - f(a));
   err = abs(c - a);
fprintf('%-6d %-14.10f %-14.10f %-14.10f %-16.10f %-16.10f\n', 0, a, b, c, f(c), err);
 for i=1:n
        if f(a)*f(c) < 0
            b=c;
        else
            a=c;
        end
        c_new = (a*f(b) - b*f(a)) / (f(b) - f(a));
        err = abs(c_new-c);
        c = c_new;
        fprintf('%-6d %-14.10f %-14.10f %-14.10f %-16.10f %-16.10f\n', i, a, b, c, f(c), err);
       if err<e
  fprintf('Root found at x = %.10f\n', c);
          break;
      end    
  end
   else
    disp('No root in this interval.')

end
