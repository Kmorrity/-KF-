clc;clear;close all ;
dbstop if error;
simTime = 500;
mtclTime = 1000;
T = 1;
F = [1 0 T 0;0 1 0 T;0 0 1 0;0 0 0 1];
q = 0.1;
Q = q*[ T^3/3, 0, T^2/2, 0; 0, T^3/3, 0 T^2/2; T^2/2, 0, T, 0; 0, T^2/2, 0, T ];
H = [1 0 0 0; 0 1 0 0];
R = diag([0.7 0.7]);

%% generate sensor measurements information


%{
figure;
hold on;
grid on;
plot( Xtrue(1, :), Xtrue(2, :), '-r.' )
Legend('Target')
%}

%{
figure;
hold on;
grid on;
plot( Xtrue(1, :), Xtrue(2, :),'-r.' )
plot( Z(1, :), Z(2, :), '-b+' );
%Legend('Target', 'Measurements')
%}

%% ini

RMSE_Rec = zeros(1, simTime, mtclTime);
RMSE_Z_Rec = zeros(1, simTime, mtclTime);

for mtclInd = 1:mtclTime
    %% generate
    Xtrue = zeros(4, simTime);
    Xtrue(:, 1) = [ 0, 0, 1, 1.5 ]';
    for simScan = 2:simTime
        Xtrue(:, simScan) = F * Xtrue(:, simScan - 1);
    end
    Xtrue = Xtrue + mvnrnd(zeros(4, 1), Q, simTime)';
    
    Z = zeros(2, simTime);
    for simScan = 1:simTime
        Z(:, simScan) = H * Xtrue(:, simScan);
    end
    Z = Z + mvnrnd(zeros(2,1), R, simTime)';
    
    %% KF
    Xupd = zeros(4, simTime);
    Xupd(:, 1) = [Z(:, 1);1;1.5];     % X后验估计 初始值
    Pupd = [R, zeros(2,2); zeros(2,2), 40^2/4*eye(2)];      % P误差协方差 初始值
    
    for simScan = 2:simScan
        Xpre = F*Xupd(:, simScan-1);
        Ppre = F*Pupd*F'+Q;
        K = Ppre*H'/(H*Ppre*H'+R);
        Xupd(:, simScan) = Xpre+K*(Z(:, simScan)-H*Xpre);
        Pupd = (eye(4)-K*H)*Ppre;

        RMSE_Z_Rec(1, simScan, mtclInd) = sum((Z(1:2, simScan)-Xtrue(1:2, simScan)).^2, 1);
        RMSE_Rec(1, simScan, mtclInd) = sum((Xupd(1:2, simScan)-Xtrue(1:2, simScan)).^2, 1);
    end
end
RMSE_Pos = mean(RMSE_Rec, 3).^0.5;
RMSE_Z_Pos = mean(RMSE_Z_Rec, 3).^0.5;

figure;
hold on;
grid on;
plot( Xupd(1, :), Xupd(2, :),'-r' )
plot( Xtrue(1, :), Xtrue(2, :),'-k' )
plot( Z(1, :), Z(2, :), '-b' );
legend('Xupd', 'Xtrue', 'Z')
xlabel('X[M]')
ylabel('Y[M]')
title('Kalman Filter')

figure;
hold on;
grid on;
plot( Xupd(3, :), '-r' )
plot( Xtrue(3, :), '-k' )
legend('Xupd_Speed', 'Xtrue_Speed')
xlabel('time[s]')
ylabel('RMSE[M]')
title('Kalman Filter Speed')

figure;
hold on;
grid on;
plot( RMSE_Pos, '-r' )
plot( RMSE_Z_Pos, '-b' );
legend( 'Est', 'Obs' )
xlabel('time[s]')
ylabel('RMSE[M]')
title('Average Position RMSE')

fprintf("Ending\n");








