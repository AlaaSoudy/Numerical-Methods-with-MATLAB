% Weighted_Gauss_Quadrature_Formula
% Student: Alaa Soudy Ibrahim Ibrahim
% Section: 5

function I = Weighted_Gauss_Quadrature_Formula()
    syms x;
    disp('Gauss Quadrature Calculator');
    disp('This program computes integral of f(x) from a to b using n=2 to 6');
    f_str = input('Enter f(x): ', 's');
    f_expr = str2sym(f_str);
    f = matlabFunction(f_expr);
    n = input('Enter n (2-6): ');
    if n < 2 || n > 6
        error('n must be between 2 and 6');
    end
    a = input('Enter a: ');
    b = input('Enter b: ');
    if b <= a
        error('b must be greater than a');
    end
    switch n
        case 2
            xi = [-0.5773502692, 0.5773502692]; wi = [1, 1];
        case 3
            xi = [-0.7745966692, 0, 0.7745966692]; wi = [0.5555555556, 0.8888888889, 0.5555555556];
        case 4
            xi = [-0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116]; wi = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451];
        case 5
            xi = [-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459]; wi = [0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851];
        case 6
            xi = [-0.9324695142, -0.6612093865, -0.2386191861, 0.2386191861, 0.6612093865, 0.9324695142]; wi = [0.1713244924, 0.3607615730, 0.4679139346, 0.4679139346, 0.3607615730, 0.1713244924];
    end
    t = ((b-a)/2)*xi + (a+b)/2;
    I = sum(wi .* f(t)) * (b-a)/2;
    try
        exact_val = integral(f, a, b, 'AbsTol', 1e-12);
        fprintf('\n∫ %s dx from %.4f to %.4f\n', f_str, a, b);
        fprintf('n = %d\n', n);
        fprintf('Gauss Result = %.10f\n', I);
        fprintf('Exact Value = %.10f\n', exact_val);
        fprintf('Error = %.2e\n', abs(I - exact_val));
    catch
        fprintf('\n∫ %s dx from %.4f to %.4f\n', f_str, a, b);
        fprintf('n = %d\n', n);
        fprintf('I = %.10f\n', I);
    end
end