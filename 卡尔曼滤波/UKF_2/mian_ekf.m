%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %程序功能：	使用扩展卡尔曼滤波器(EKF)估计平抛物体的运动
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clc;
 clear;
 close all;
kx = .01; ky = .01; 	% 阻尼系数
    g = 10; 		% 重力
    t = 10; 		% 仿真时间
    Ts = 0.1; 		% 采样周期
    len = fix(t/Ts);    % 仿真步数
    m = 10;         %质量
    %（真实轨迹模拟）
    dax = 0.1; day = 0.1;  % 系统噪声
    X = zeros(len,4); X(1,:) = [0, 50, 500, 0]; % 状态模拟的初值
    for k=2:len
        x = X(k-1,1); vx = X(k-1,2); y = X(k-1,3); vy = X(k-1,4); 
        x = x + vx*Ts;
        vx = vx + (-kx*vx^2)*Ts+dax*randn(1);
        y = y + vy*Ts;
        vy = vy + ((ky*vy^2-m*g)/m)*Ts+day*randn(1);
        X(k,:) = [x, vx, y, vy];
    end
    %figure(1), hold off, plot(X(:,2),'-b'), grid on
    
    figure(1), hold off, plot(X(:,1),X(:,3),'-b'), grid on   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 构造量测量
    mrad = 0.001;
    dr = 10; dafa = 10*mrad; % 量测噪声
    for k=1:len
        r = sqrt(X(k,1)^2+X(k,3)^2) + dr*randn(1,1);
        a = atan(X(k,1)/X(k,3)) + dafa*randn(1,1);
        Z(k,:) = [r, a];
    end
    figure(1), hold on, plot(Z(:,1).*sin(Z(:,2)), Z(:,1).*cos(Z(:,2)),'+k')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ekf 滤波
    Qk = diag([0; dax; 0; day])^2;
    Rk = diag([dr; dafa])^2;
    Xk = zeros(4,1);
    Pk = 100*eye(4);
    X_est = X;
    for k=1:len
        Ft = JacobianF(X(k,:), kx, ky, g);
        Hk = JacobianH(X(k,:));
        fX = fff(X(k,:), kx, ky, g, Ts);
        hfX = hhh(fX, Ts);
        [Xk, Pk, Kk] = ekf(eye(4)+Ft*Ts, Qk, fX, Pk, Hk, Rk, Z(k,:)'-hfX);
        X_est(k,:) = Xk';
    end
    figure(1), plot(X_est(:,1),X_est(:,3), 'r+')
    xlabel('X'); ylabel('Y'); title('ekf simulation');
    legend('real', 'measurement', 'ekf estimated');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
