function [pos,vel,posPlot,quatPlot,gyrX,gyrY,gyrZ,accX,accY,accZ,acc,time,acc_magFilt,stationary,samplePeriod,ang]=script2(str,sensor,filt,filtrobajo)
cd ../
addpath('Quaternions');
addpath('ximu_matlab_library');
[A,~,~,t,sensores,G,~,~]=cadera(str);
l=1;
if sensor>sensores
    sensor=sensores;
elseif sensor<1
    sensor=1;
end

for k=1:length(t)
    if A(1,k,sensor)==0
        l=l+1;
    else
        break
    end
end

m=1;
for k=1:length(t)
    if A(1,length(t)-k+1,sensor)==0
        m=m+1;
    else
        break
    end
end
m=m+1;
samplePeriod=0.00675;
filtro=filt;
time=t(l:length(t)-m);
accX=A(1,l:length(t)-m,sensor)';
accY=A(2,l:length(t)-m,sensor)';
accZ=A(3,l:length(t)-m,sensor)';
gyrX=G(1,l:length(t)-m,sensor)';
gyrY=G(2,l:length(t)-m,sensor)';
gyrZ=G(3,l:length(t)-m,sensor)';

% Manually frame data

% startTime = 0;
% stopTime = 10;
% -------------------------------------------------------------------------
% Detect stationary periods

% Compute accelerometer magnitude
acc_mag = sqrt(accX.*accX + accY.*accY + accZ.*accZ);
% [f,espectro]=fftt(acc_mag,time);
% plot(f,espectro)
% A=espectro>0.1;
% f=f(A)
% pause
% HP filter accelerometer data
filtCutOff = filtrobajo;
[b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'high');
acc_magFilt = filtfilt(b, a, acc_mag);

% Compute absolute value
acc_magFilt = abs(acc_magFilt);

% LP filter accelerometer data
filtCutOff = 3;
[b, a] = butter(1, (2*filtCutOff)/(1/samplePeriod), 'low');
acc_magFilt = filtfilt(b, a, acc_magFilt);

% Threshold detection
stationary = acc_magFilt < filtro;

% -------------------------------------------------------------------------
% Plot data raw sensor data and stationary periods

% figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Sensor Data');
% ax(1) = subplot(2,1,1);
%     hold on;
%     plot(time, gyrX, 'r');
%     plot(time, gyrY, 'g');
%     plot(time, gyrZ, 'b');
%     title('Gyroscope');
%     xlabel('Time (s)');
%     ylabel('Angular velocity (^\circ/s)');
%     legend('X', 'Y', 'Z');
%     hold off;
% ax(2) = subplot(2,1,2);
%     hold on;
%     plot(time, accX, 'r');
%     plot(time, accY, 'g');
%     plot(time, accZ, 'b');
%     plot(time, acc_magFilt, ':k');
%     plot(time, stationary, 'k', 'LineWidth', 2);
%     title('Accelerometer');
%     xlabel('Time (s)');
%     ylabel('Acceleration (g)');
%     legend('X', 'Y', 'Z', 'Filtered', 'Stationary');
%     hold off;
% linkaxes(ax,'x');

% -------------------------------------------------------------------------
% Compute orientation

quat = zeros(length(time), 4);
AHRSalgorithm = AHRS('SamplePeriod', 1/256, 'Kp', 1, 'KpInit', 1);

% Initial convergence
initPeriod = 2;
indexSel = 6 : find(sign(time-(time(1)+initPeriod))+1, 1);
[li ai]=size(indexSel);
for i = 1:2000
    AHRSalgorithm.UpdateIMU([0 0 0], [sum(accX(indexSel))/ai sum(accY(indexSel))/ai sum(accZ(indexSel))/ai]);
end

% For all data
for k = 1:length(time)
    if(stationary(k))
        AHRSalgorithm.Kp = 0.5;
    else
        AHRSalgorithm.Kp = 0;
    end
    AHRSalgorithm.UpdateIMU(deg2rad([gyrX(k) gyrY(k) gyrZ(k)]), [accX(k) accY(k) accZ(k)]);
    quat(k,:) = AHRSalgorithm.Quaternion;
end

% -------------------------------------------------------------------------
% Compute translational accelerations

% Rotate body accelerations to Earth frame
acc = quaternRotate([accX accY accZ], quaternConj(quat));

% % Remove gravity from measurements
% acc = acc - [zeros(length(time), 2) ones(length(time), 1)];     % unnecessary due to velocity integral drift compensation

% Convert acceleration measurements to m/s/s
acc = acc * 9.81;

% Plot translational accelerations
% figure('Position', [9 39 900 300], 'NumberTitle', 'off', 'Name', 'Accelerations');
% hold on;
% plot(time, acc(:,1), 'r');
% plot(time, acc(:,2), 'g');
% plot(time, acc(:,3), 'b');
% title('Acceleration');
% xlabel('Time (s)');
% ylabel('Acceleration (m/s/s)');
% legend('X', 'Y', 'Z');
% hold off;

% -------------------------------------------------------------------------
% Compute translational velocities

acc(:,3) = acc(:,3) - 9.81;

% Integrate acceleration to yield velocity
vel = zeros(size(acc));
for k = 2:length(vel)
    vel(k,:) = vel(k-1,:) + acc(k,:) * samplePeriod;
    if(stationary(k) == 1)
        vel(k,:) = [0 0 0];     % force zero velocity when foot stationary
    end
end


% Compute integral drift during non-stationary periods
velDrift = zeros(size(vel));
stationaryStart = find([0; diff(stationary)] == -1);
stationaryEnd = find([0; diff(stationary)] == 1);
for i = 1:numel(stationaryEnd)-1
    driftRate = vel(stationaryEnd(i)-1, :) / (stationaryEnd(i) - stationaryStart(i));
    enum = 1:(stationaryEnd(i) - stationaryStart(i));
    drift = [enum'*driftRate(1) enum'*driftRate(2) enum'*driftRate(3)];
    velDrift(stationaryStart(i):stationaryEnd(i)-1, :) = drift;
end

% Remove integral drift
vel = vel - velDrift;

% Plot translational velocity
% figure('Position', [9 39 900 300], 'NumberTitle', 'off', 'Name', 'Velocity');
% hold on;
% plot(time, vel(:,1), 'r');
% plot(time, vel(:,2), 'g');
% plot(time, vel(:,3), 'b');
% title('Velocity');
% xlabel('Time (s)');
% ylabel('Velocity (m/s)');
% legend('X', 'Y', 'Z');
% hold off;

% -------------------------------------------------------------------------
% Compute translational position

% Integrate velocity to yield position
pos = zeros(size(vel));
for k = 2:length(pos)
    pos(k,:) = pos(k-1,:) + vel(k,:) * samplePeriod;    % integrate velocity to yield position
end

gyr=[gyrX gyrY gyrZ];
ang= zeros(size(gyr));
for k = 2:length(vel)
    ang(k,:) = ang(k-1,:) + gyr(k,:) * samplePeriod;
    if(stationary(k) == 1)
        ang(k,:) =[0 0 0];     % force zero velocity when foot stationary
    end
end
% % Plot translational position
% figure('Position', [9 39 900 600], 'NumberTitle', 'off', 'Name', 'Position');
% hold on;
% plot(time, pos(:,1), 'r');
% plot(time, pos(:,2), 'g');
% plot(time, pos(:,3), 'b');
% title('Position');
% xlabel('Time (s)');
% ylabel('Position (m)');
% legend('X', 'Y', 'Z');
% hold off;
% 
% % -------------------------------------------------------------------------
% % Plot 3D foot trajectory
% 
% % % Remove stationary periods from data to plot
% posPlot = pos(find(~stationary), :);
% quatPlot = quat(find(~stationary), :);
posPlot = pos;
quatPlot = quat;

% Extend final sample to delay end of animation
extraTime = 20;
onesVector = ones(floor(extraTime*(1/samplePeriod))+1, 1);
posPlot = [posPlot; [posPlot(end, 1)*onesVector, posPlot(end, 2)*onesVector, posPlot(end, 3)*onesVector]];
quatPlot = [quatPlot; [quatPlot(end, 1)*onesVector, quatPlot(end, 2)*onesVector, quatPlot(end, 3)*onesVector, quatPlot(end, 4)*onesVector]];
% 
% % Create 6 DOF animation
% SamplePlotFreq = 4;
% Spin = 120;
% SixDofAnimation(posPlot, quatern2rotMat(quatPlot), ...
%                 'SamplePlotFreq', SamplePlotFreq, 'Trail', 'All', ...
%                 'Position', [9 39 1280 768], 'View', [(100:(Spin/(length(posPlot)-1)):(100+Spin))', 10*ones(length(posPlot), 1)], ...
%                 'AxisLength', 0.1, 'ShowArrowHead', false, ...
%                 'Xlabel', 'X (m)', 'Ylabel', 'Y (m)', 'Zlabel', 'Z (m)', 'ShowLegend', false, ...
%                 'CreateAVI', false, 'AVIfileNameEnum', false, 'AVIfps', ((1/samplePeriod) / SamplePlotFreq));
 end