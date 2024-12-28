clearvars; 
clc; 
close all;
%Loading the true data belonging to the target 
load("multi_target_data.mat");

T = 1;
t = 0:T:100;

xrange = 10e2;
yrange = 10e2;
V = xrange * yrange;

measurement_std = 10;
measurement_mu = [0;0];

beta_FA = 1e-5; %clutter intensity

P_d = 0.9;
P_s = 0.99;

A = [eye(2),T*eye(2);0*eye(2),eye(2)];
G = [T^2/2*eye(2);T*eye(2)];
C = [eye(2),zeros(2)];
Q = [2^2,0;0,2^2];
R = measurement_std^2*eye(2);

x1data = data{1,1};
x2data = data{2,1};
x3data = data{3,1};
x4data = data{4,1};

trueTargets = cell(3,4);
trueTargets{1,1} = x1data(2,:);
trueTargets{2,1} = x1data(3,:);
trueTargets{3,1} = x1data(1,:);
trueTargets{1,2} = x2data(2,:);
trueTargets{2,2} = x2data(3,:);
trueTargets{3,2} = x2data(1,:);
trueTargets{1,3} = x3data(2,:);
trueTargets{2,3} = x3data(3,:);
trueTargets{3,3} = x3data(1,:);
trueTargets{1,4} = x4data(2,:);
trueTargets{2,4} = x4data(3,:);
trueTargets{3,4} = x4data(1,:);



figure;
plot(x1data(2,:),x1data(3,:),"b.");
hold on;
plot(x2data(2,:),x2data(3,:),"g.");
plot(x3data(2,:),x3data(3,:),"r.");
plot(x4data(2,:),x4data(3,:),"y.");
xlim([0,xrange]);
ylim([0,yrange]);
xlabel("x");
ylabel("y");
title("Standard Illustration: Missing Time Data");
grid on;

figure;
sgtitle("MTT Illustration");
subplot(2,1,1);
plot(x1data(1,:),x1data(2,:),"b.");
hold on;
plot(x2data(1,:),x2data(2,:),"g.");
plot(x3data(1,:),x3data(2,:),"r.");
plot(x4data(1,:),x4data(2,:),"y.");
xlim([0,300]);
ylim([0,xrange]);
xlabel("Time (s)");
ylabel("x");
grid on;

subplot(2,1,2);
plot(x1data(1,:),x1data(3,:),"b.");
hold on;
plot(x2data(1,:),x2data(3,:),"g.");
plot(x3data(1,:),x3data(3,:),"r.");
plot(x4data(1,:),x4data(3,:),"y.");
xlim([0,300]);
ylim([0,yrange]);
xlabel("Time (s)");
ylabel("y");
grid on;


% clutter generation

false_alarms = cell(3,length(t));
for i= 1:length(t)
    m_k = poissrnd(beta_FA*V);
    FA_xvalues = xrange*rand(1,m_k);
    FA_yvalues = yrange*rand(1,m_k);
    false_alarms{1,i} = t(i);
    false_alarms{2,i} = FA_xvalues;
    false_alarms{3,i} = FA_yvalues;
end

measurements1 = zeros(3,length(t));
measurements2 = zeros(3,length(t));
measurements3 = zeros(3,length(t));
measurements4 = zeros(3,length(t));

for i = 1:length(t)
    u1 = rand();
    u2 = rand();
    u3 = rand();
    u4 = rand();
    k = find(x1data(1,:) == t(i),1);
    if ~isempty(k)
        if u1<=P_d
        % detection
        v_k = mvnrnd(measurement_mu,R)';
        measurement = x1data(2:3,k) +v_k;
        measurements1(1,i) = t(i);
        measurements1(2:3,i) = measurement;
        else
        % miss
        measurements1(1,i) = t(i);
        measurements1(2,i) = inf;
        measurements1(3,i) = inf;
        end
    else
        % target has not yet begun its track
        measurements1(1,i) = t(i);
        measurements1(2,i) = inf;
        measurements1(3,i) = inf;
    end

    k = find(x2data(1,:) == t(i),1);
    if ~isempty(k)
        if u2<=P_d
        % detection
        v_k = mvnrnd(measurement_mu,R)';
        measurement = x2data(2:3,k) +v_k;
        measurements2(1,i) = t(i);
        measurements2(2:3,i) = measurement;
        else
        % miss
        measurements2(1,i) = t(i);
        measurements2(2,i) = inf;
        measurements2(3,i) = inf;
        end
    else
        % target has not yet begun its track
        measurements2(1,i) = t(i);
        measurements2(2,i) = inf;
        measurements2(3,i) = inf;
    end

    k = find(x3data(1,:) == t(i),1);
    if ~isempty(k)
        if u3<=P_d
        % detection
        v_k = mvnrnd(measurement_mu,R)';
        measurement = x3data(2:3,k) +v_k;
        measurements3(1,i) = t(i);
        measurements3(2:3,i) = measurement;
        else
        % miss
        measurements3(1,i) = t(i);
        measurements3(2,i) = inf;
        measurements3(3,i) = inf;
        end
    else
        % target not yet begun its track
        measurements3(1,i) = t(i);
        measurements3(2,i) = inf;
        measurements3(3,i) = inf;
    end
    k = find(x4data(1,:) == t(i),1);
    if ~isempty(k)
        if u4<=P_d
        % detection
        v_k = mvnrnd(measurement_mu,R)';
        measurement = x4data(2:3,k) +v_k;
        measurements4(1,i) = t(i);
        measurements4(2:3,i) = measurement;
        else
        % miss
        measurements4(1,i) = t(i);
        measurements4(2,i) = inf;
        measurements4(3,i) = inf;
        end
    else
        % target not yet begun its track
        measurements4(1,i) = t(i);
        measurements4(2,i) = inf;
        measurements4(3,i) = inf;
    end
end

detections = cell(3,length(t));

for i = 1:length(t)
    detections{1,i} = t(i);
    detections{2,i} = false_alarms{2,i};
    detections{3,i} = false_alarms{3,i};
    if measurements1(2,i) ~= inf
        detections{2,i} = [detections{2,i} measurements1(2,i)];
        detections{3,i} = [detections{3,i} measurements1(3,i)];
    end
    if measurements2(2,i) ~= inf
        detections{2,i} = [detections{2,i} measurements2(2,i)];
        detections{3,i} = [detections{3,i} measurements2(3,i)];
    end
    if measurements3(2,i) ~= inf
        detections{2,i} = [detections{2,i} measurements3(2,i)];
        detections{3,i} = [detections{3,i} measurements3(3,i)];
    end
    if measurements4(2,i) ~= inf
        detections{2,i} = [detections{2,i} measurements4(2,i)];
        detections{3,i} = [detections{3,i} measurements4(3,i)];
    end
end

estimated_states = GM_PHD_Filter(detections,t,P_d,P_s,beta_FA,A,G,C,Q,R,V);

figure;

for k = t
    tracks = estimated_states{k+1};
    for i = 1:length(tracks)
        current_track = tracks{i};
        scatter(current_track(1),current_track(2),'blue');
        hold on;
    end
end

xlim([0,xrange]);
ylim([0,yrange]);
xlabel("x");
ylabel("y");
title("Estimated States");
grid on;


figure;

for k = t
    xvalues = detections{2,k+1};
    yvalues = detections{3,k+1};

    for i = 1:length(xvalues)
        scatter(xvalues(i),yvalues(i),'black');
        hold on;
    end
end

xlim([0,xrange]);
ylim([0,yrange]);
xlabel("x");
ylabel("y");
title("Detections");
grid on;


figure;
hold on;
for k = t
    xvalues = detections{2,k+1};
    yvalues = detections{3,k+1};

    for i = 1:length(xvalues)
        scatter(xvalues(i),yvalues(i),'black');
    end

    tracks = estimated_states{k+1};
    for i = 1:length(tracks)
        current_track = tracks{i};
        scatter(current_track(1),current_track(2),'blue');
        hold on;
    end

    drawnow;
    pause(0.1);

end


