clc;
clear;
close all;
%warning off;

%{
X_{k+1} = F_k * X_k + V_k;
Z_{k+1} = h(X_{k+1}, W_K);
Z_{k+1} = [r_{k+1}, \theta_{k+1}];
%} 

%% set params
TarDims = 4; % longitudianal and lateral position and velocity
MeaDims = 2; % complete linear measurement - range and theta
T = 1; % simple interval
simTime = 1e2; % simulation total time
mtclTime = 5e2; % simulation total mtcl times
F = kron([1 T; 0, 1], eye(2)); % target state transform matrix
% H = eye(4); % sensor observe matrix
% H = [];
q = 1e-2;
Q = q * kron([T^3/3, T^2/2; T^2/2, T], eye(2)); % process noise
R = diag([10, deg2rad(3)] .^ 2); % measurements noise
X0 = [1e3, 2e2, 20, -10]; % target initialization state
Vmax = 40; % target maximum linear velocity

%% generate target trajectory and sensor measurements
XTrue = zeros(TarDims, simTime);
ZMeas = zeros(MeaDims, simTime);
XTrue(:, 1) = X0;
ZMeas(:, 1) = [norm(X0(1 : 2)); atan2(X0(2), X0(1))];
%% generate target moving path and observations
for simScan = 2 : T : simTime
    XTrue(:, simScan) = F * XTrue(:, simScan - 1);
    ZMeas(:, simScan) = [norm(XTrue(1 : 2, simScan)), atan2(XTrue(2, simScan), XTrue(1, simScan))];
end
% add process and measurement noise
XTrue = repmat(XTrue, [1, 1, mtclTime]);
ZMeas = repmat(ZMeas, [1, 1, mtclTime]);
for mtclInd = 1 : mtclTime
    XTrue(:, :, mtclInd) = XTrue(:, :, mtclInd) + mvnrnd(zeros(TarDims, 1), Q, simTime)';
    ZMeas(:, :, mtclInd) = ZMeas(:, :, mtclInd) + mvnrnd(zeros(MeaDims, 1), R, simTime)';
end
%{
% plot figures to identify target trajectory and measurements
Z_XY = [ZMeas(1, :, 1) .* cos(ZMeas(2, :, 1)); ZMeas(1, :, 1) .* sin(ZMeas(2, :, 1))];
figure
hold on;
grid on;
plot(XTrue(1, :, 1), XTrue(2, :, 1), '-r.');
plot(Z_XY(1, :, 1), Z_XY(2, :, 1), 'b*', 'markersize', 10);

legend('Taget Trajectory', 'Sensor Measurements', 'Location', 'best');
set(get(gca, 'XLabel'), 'String', 'Xaxis-[m]');
set(get(gca, 'YLabel'), 'String', 'Yaxis-[m]');
saveas(gcf, './Trajectory and measurements', 'fig');
saveas(gcf, './Trajectory and measurements', 'bmp');
%}


%% filter
Xupd = zeros(TarDims, simTime, mtclTime);
Pupd = zeros(TarDims, TarDims, simTime, mtclTime);
ZCon = zeros(2, simTime, mtclTime);
for mtclInd = 1 : mtclTime
    % initialize estimate mean and covariance
    % Xupd(:, 1, mtclInd) = ZMeas(:, 1, mtclInd);
    Xupd(:, 1, mtclInd) = [ ZMeas(1, 1, mtclInd) * cos(ZMeas(2, 1, mtclInd));
                            ZMeas(1, 1, mtclInd) * sin(ZMeas(2, 1, mtclInd));
                            0; 0];
    Pupd(:, :, 1, mtclInd) = 10 ^ 2 * eye(4);    % 无偏转换
    
    ZCon(1, 1, mtclInd) = ZMeas(1, 1, mtclInd) * cos(ZMeas(2, 1, mtclInd));
    ZCon(2, 1, mtclInd) = ZMeas(1, 1, mtclInd) * sin(ZMeas(2, 1, mtclInd));
    for simScan = 2 : simTime
        H = zeros(MeaDims, TarDims);
        Xpre = F * Xupd(:, simScan - 1, mtclInd);
        Ppre = F * Pupd(:, :, simScan - 1, mtclInd) * F' + Q;
        % cal H
        range_pre = norm(Xpre(1 : 2));
        x_pre = Xpre(1);
        y_pre = Xpre(2);
        H(1, 1) = x_pre / range_pre;
        H(1, 2) = y_pre / range_pre;
        H(2, 1) = -1 * y_pre / (range_pre ^ 2);
        H(2, 2) = x_pre / (range_pre ^ 2);

        % rest
        Zpre = [norm(Xpre(1 : 2)); atan2(Xpre(2), Xpre(1))];
        Sk = H * Ppre * H' + R;
        K = Ppre * H' / Sk;

        Xupd(:, simScan, mtclInd) = Xpre + K * (ZMeas(:, simScan, mtclInd) - Zpre);
        Pupd(:, :, simScan, mtclInd) = (eye(TarDims) - K * H) * Ppre;

        ZCon(1, simScan, mtclInd) = ZMeas(1, simScan, mtclInd) * cos(ZMeas(2, simScan, mtclInd));
        ZCon(2, simScan, mtclInd) = ZMeas(1, simScan, mtclInd) * sin(ZMeas(2, simScan, mtclInd));
    end
end

%% calculate the performance
RMSETemp = Xupd - XTrue;
% RMSEZTemp = ZMeas - XTrue;
RMSEZTemp = ZCon - XTrue(1 : 2, :, :);
RMSE = mean(RMSETemp .^ 2, 3) .^ 0.5;
RMSEZ = mean(RMSEZTemp .^ 2, 3) .^ 0.5;

%% show the performance
figure
hold on;
grid on;
plot(Xupd(1, :, 1), Xupd(2, :, 1), '-k.');
plot(ZCon(1, :, 1), ZCon(2, :, 1), 'b*');
plot(XTrue(1, :, 1), XTrue(2, :, 1), '-r.')
legend('Esti', 'Meas', 'True');
set(get(gca, 'XLabel'), 'String', 'Xaxis-[m]');
set(get(gca, 'YLabel'), 'String', 'Yaxis-[m]');
title('Target True Trajectory, Measurements and Corresponding KF estimation')
saveas(gcf, './Trajectory, measurements and estimations', 'fig');
saveas(gcf, './Trajectory, measurements and estimations', 'bmp');

for plotInd = 1 : 2
    figure
    hold on;
    grid on;
    plot(RMSE(plotInd, :), '-r.')
    plot(RMSEZ(plotInd, :), '-k.')
    legend('Esti', 'Meas')
    set(get(gca, 'XLabel'), 'String', 'time-[s]');
    switch plotInd
    case 1
        set(get(gca, 'YLabel'), 'String', 'RMSE-[m]');
        title('Position X')        
        saveas(gcf, './Lateral Position', 'fig');
        saveas(gcf, './Lateral Position', 'bmp');
    case 2
        set(get(gca, 'YLabel'), 'String', 'RMSE-[m]');
        title('Position Y')
        saveas(gcf, './Longitudinal Position', 'fig');
        saveas(gcf, './Longitudinal Position', 'bmp');
    end
end


for plotInd = 1 : 2
    figure
    grid on;
    plot(RMSE(plotInd + 2, :), '-r.')
    legend('Esti')
    set(get(gca, 'XLabel'), 'String', 'time-[s]');
    set(get(gca, 'YLabel'), 'String', 'RMSE-[m/s]');
    
    switch plotInd
    case 1
        title('Velocity X')        
        saveas(gcf, './Lateral Velocity', 'fig');
        saveas(gcf, './Lateral Velocity', 'bmp');
    case 2
        title('Velocity Y')
        saveas(gcf, './Longitudinal Velocity', 'fig');
        saveas(gcf, './Longitudinal Velocity', 'bmp');
    end

end



%{

H = \frac{Z_k}{X_k} = zeros(2, 4);
H(1, 1) = \frac{range_k}{x_k} = sqrt{x^2 + y^2}{x} = x/r;
H(1, 2) = y / r;
H(1, 3) = 0;
H(1, 4) = 0;
H(2, 1) = -1 * y / (r ^ 2);
H(2, 2) = x / r ^ 2;
H(2, 3) = 0;
H(2, 4) = 0;

%}