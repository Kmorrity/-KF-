function F = JacobianF(X, kx, ky, g) % ϵͳ״̬�ſɱȺ���
    vx = X(2); vy = X(4); 
    F = zeros(4,4);
    F(1,2) = 1;
    F(2,2) = -2*kx*vx;
    F(3,4) = 1;
    F(4,4) = 2*ky*vy;
