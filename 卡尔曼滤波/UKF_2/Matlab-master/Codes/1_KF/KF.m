clc;clear;close all;

dbstop if error;
warning off;

TarDims = 4; % longitudianal and lateral position and velocity
MeaDims = 4; % complete linear measurement
T = 1; % simple interval
simTime = 1e2; % simulation total time
mtclTime = 5e2; % simulation total mtcl times
F = kron([1 T; 0, 1], eye(2)); % target state transform matrix
H = eye(4); % sensor observe matrix
q = 1e-2;
Q = q * kron([T^3/3, T^2/2; T^2/2, T], eye(2)); % process noise
R = diag([10, 10, 5, 3] .^ 2); % measurements noise
X0 = [1e3, 2e2, 20, -10]; % target initialization state
Vmax = 40; % target maximum linear velocity

%% preallocated space for matrix
XTrue = zeros(TarDims, simTime);
ZMeas = zeros(MeaDims, simTime);
XTrue(:, 1) = X0;
ZMeas(:, 1) = H * XTrue(:, 1);
%% generate target moving path and observations
for simScan = 2 : T : simTime
    XTrue(:, simScan) = F * XTrue(:, simScan - 1);
    ZMeas(:, simScan) = H * XTrue(:, simScan);
end
% add process and measurement noise
XTrue = repmat(XTrue, [1, 1, mtclTime]);
ZMeas = repmat(ZMeas, [1, 1, mtclTime]);
for mtclInd = 1 : mtclTime
    XTrue(:, :, mtclInd) = XTrue(:, :, mtclInd) + mvnrnd(zeros(TarDims, 1), Q, simTime)';
    ZMeas(:, :, mtclInd) = ZMeas(:, :, mtclInd) + mvnrnd(zeros(MeaDims, 1), R, simTime)';
end
% %{
% plot figures to identify target trajectory and measurements
figure
hold on;
grid on;
plot(XTrue(1, :, 1), XTrue(2, :, 1), '-r.');
plot(ZMeas(1, :, 1), ZMeas(2, :, 1), 'b*', 'markersize', 10);
legend('Taget Trajectory', 'Sensor Measurements', 'Location', 'best');
set(get(gca, 'XLabel'), 'String', 'Xaxis-[m]');
set(get(gca, 'YLabel'), 'String', 'Yaxis-[m]');
saveas(gcf, './Trajectory and measurements', 'fig');
saveas(gcf, './Trajectory and measurements', 'bmp');
%}

%% filtering
Xupd = zeros(TarDims, simTime, mtclTime);
Pupd = zeros(TarDims, TarDims, simTime, mtclTime);
for mtclInd = 1 : mtclTime
    % initialize estimate mean and covariance
    Xupd(:, 1, mtclInd) = ZMeas(:, 1, mtclInd);
    Pupd(:, :, 1, mtclInd) = R;    
    for simScan = 2 : simTime
        Xpre = F * Xupd(:, simScan - 1, mtclInd);
        Ppre = F * Pupd(:, :, simScan - 1, mtclInd) * F' + Q;

        Zpre = H * Xpre;
        Sk = H * Ppre * H' + R;
        K = Ppre * H' / Sk;

        Xupd(:, simScan, mtclInd) = Xpre + K * (ZMeas(:, simScan, mtclInd) - Zpre);
        Pupd(:, :, simScan, mtclInd) = (eye(TarDims) - K * H) * Ppre;
    end
end

RMSETemp = Xupd - XTrue;
RMSEZTemp = ZMeas - XTrue;
RMSE = mean(RMSETemp .^ 2, 3) .^ 0.5;
RMSEZ = mean(RMSEZTemp .^ 2, 3) .^ 0.5;

figure
hold on;
grid on;
plot(Xupd(1, :, 1), Xupd(2, :, 1), '-k.');
plot(ZMeas(1, :, 1), ZMeas(2, :, 1), 'b*');
plot(XTrue(1, :, 1), XTrue(2, :, 1), '-r.')
legend('Esti', 'Meas', 'True');
set(get(gca, 'XLabel'), 'String', 'Xaxis-[m]');
set(get(gca, 'YLabel'), 'String', 'Yaxis-[m]');
title('Target True Trajectory, Measurements and Corresponding KF estimation')
saveas(gcf, './Trajectory, measurements and estimations', 'fig');
saveas(gcf, './Trajectory, measurements and estimations', 'bmp');

for plotInd = 1 : TarDims
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
    case 3
        set(get(gca, 'YLabel'), 'String', 'RMSE-[m/s]');
        title('Velocity X')
        saveas(gcf, './Lateral Velocity', 'fig');
        saveas(gcf, './Lateral Velocity', 'bmp');
    case 4
        set(get(gca, 'YLabel'), 'String', 'RMSE-[m/s]');
        title('Velocity Y')
        saveas(gcf, './Longitudinal Velocity', 'fig');
        saveas(gcf, './Longitudinal Velocity', 'bmp');
    end
end