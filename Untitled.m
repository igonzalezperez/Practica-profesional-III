close all
str='IvanCaminando.csv';
sensor=2;
filt=0.45;
filtrobajo=0.004;
[pos,vel,posPlot,quatPlot,gyrX,gyrY,gyrZ,accX,accY,accZ,acc,time,acc_magFilt,stationary,samplePeriod]=script2(str,sensor,filt,filtrobajo);
hold on
gyr=[gyrX gyrY gyrZ];
ang= zeros(size(gyr));
for k = 2:length(vel)
    ang(k,:) = ang(k-1,:) + gyr(k,:) * samplePeriod;
    if(stationary(k) == 1)
        ang(k,:) =[0 0 0];     % force zero velocity when foot stationary
    end
end