function fX = fff(X, kx, ky, g, Ts) % 系统状态非线性函数
    x = X(1); vx = X(2); y = X(3); vy = X(4); 
    x1 = x + vx*Ts;
    vx1 = vx + (-kx*vx^2)*Ts;
    y1 = y + vy*Ts;
    vy1 = vy + (ky*vy^2-g)*Ts;
    fX = [x1; vx1; y1; vy1];
