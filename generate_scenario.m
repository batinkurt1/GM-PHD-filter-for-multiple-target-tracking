T = 1; 
timesteps = 0:T:100; 
xrange = 10e2;
yrange = 10e2;


A = [eye(2), T*eye(2); zeros(2), eye(2)];
G = [T^2/2*eye(2); T*eye(2)];
C = [1, 0, 0, 0; 0, 1, 0, 0];
Q = 2^2 * eye(2);



data = cell(5, 1);


birth_time1 = 5;
death_time1 = 77;
initial_state1 = [500; 500; 10*(rand-0.5); 10*(rand-0.5)]; % Initial state: [x; y; vx; vy]
data{1} = simulate_target(initial_state1, A, G, Q, birth_time1, death_time1);


birth_time2 = 45;
death_time2 = 63;
initial_state2 = data{1}(2:5, birth_time2-birth_time1+1); % Initial state of target 2
initial_state2(3:4) = [initial_state2(4); -initial_state2(3)]; % Perpendicular direction
data{2} = simulate_target(initial_state2, A, G, Q, birth_time2, death_time2);


birth_time3 = 10;
death_time3 = 86;
initial_state3 = [250; 400; 10*(rand-0.5); 10*(rand-0.5)];
data{3} = simulate_target(initial_state3, A, G, Q, birth_time3, death_time3);


birth_time4 = 56;
death_time4 = 90;
initial_state4 = data{3}(2:5, birth_time4-birth_time3+1);
initial_state4(3:4) = [initial_state4(4); -initial_state4(3)]; % Perpendicular direction
data{4} = simulate_target(initial_state4, A, G, Q, birth_time4, death_time4);




save('multi_target_data.mat', 'data');


figure; hold on;
colors = ['r', 'g', 'b', 'm','c'];
for i = 1:4
    plot(data{i}(2, :), data{i}(3, :), 'Color', colors(i), 'DisplayName', ['Target ' num2str(i)]);
end
legend show;
xlabel('X position');
ylabel('Y position');
title('Multiple Target Scenario');
axis([0 xrange 0 yrange]);
grid on;


function target_data = simulate_target(initial_state, A, G, Q, birth_time, death_time)
    timesteps = birth_time:death_time;
    state = initial_state;
    target_data = zeros(5, length(timesteps));
    for t = 1:length(timesteps)
        state = A * state + G * mvnrnd([0; 0], Q)';
        %state = A * state;
        target_data(:, t) = [timesteps(t); state];
    end
end
