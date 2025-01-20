% MATLAB Code to Visualize Lid-Driven Cavity Results with Contour and Quiver

% Load the results from the CSV file
data = readmatrix('2D_results_square.csv');

% Extract columns from the data
x = data(:, 1); % x-coordinates
y = data(:, 2); % y-coordinates
u = data(:, 3); % u-velocity
v = data(:, 4); % v-velocity

% Reshape the data into grid form
N = sqrt(length(x)) - 1; % Assuming a square grid
x = reshape(x, N+1, N+1);
y = reshape(y, N+1, N+1);
u = reshape(u, N+1, N+1);
v = reshape(v, N+1, N+1);

% Compute speed magnitude
speed = sqrt(u.^2 + v.^2);

% Plot the speed magnitude (contour plot) with velocity vectors
figure;
contourf(x, y, speed, 20, 'LineStyle', 'none'); % Contour plot of speed
colorbar;
colormap("jet");
hold on; % Overlay the quiver plot
quiver(x, y, u, v, 'k'); % Velocity vectors (black arrows)
title('Speed Magnitude Contour with Velocity Vectors');
xlabel('X');
ylabel('Y');
axis equal;
grid on;
hold off;

% Find all rows where y = 0.5 for v-velocity plot
y_target = 0.5;
tolerance = 1e-6; % Numerical tolerance for floating-point comparison
indices_y = find(abs(y - y_target) < tolerance); % Rows where y = 0.5
x_at_y_half = x(indices_y);
v_at_y_half = v(indices_y);

% Plot v-velocity along y = H/2
figure;
plot(x_at_y_half, v_at_y_half, 'b-', 'LineWidth', 2); % Simulated data
hold on;
title('Re=1000, V-velocity at y=H/2');
xlabel('x/L');
ylabel('V/V0');
legend('Simulated', 'Location', 'Best');
grid on;

% Find all rows where x = 0.5 for u-velocity plot
x_target = 0.5;
indices_x = find(abs(x - x_target) < tolerance); % Rows where x = 0.5
y_at_x_half = y(indices_x);
u_at_x_half = u(indices_x);

% Plot u-velocity along x = H/2
figure;
plot(u_at_x_half, y_at_x_half, 'r-', 'LineWidth', 2); % Simulated data
hold on;
title('Re=1000, U-velocity at x=H/2');
xlabel('U/U0');
ylabel('y/L');
legend('Simulated', 'Location', 'Best');
grid on;

