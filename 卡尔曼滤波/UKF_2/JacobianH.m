function H = JacobianH(X) % 量测雅可比函数
    x = X(1); y = X(3);
    H = zeros(2,4);
    r = sqrt(x^2+y^2);
    H(1,1) = 1/r; H(1,3) = 1/r;
    xy2 = 1+(x/y)^2;
    H(2,1) = 1/xy2*1/y; H(2,3) = 1/xy2*x*(-1/y^2);
