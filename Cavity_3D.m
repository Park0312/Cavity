% Load the CSV file containing x, y, z, u, v, w data
data = readmatrix('3D_results_100seconds.csv');

% Extract columns for coordinates and velocities
x = data(:, 1);
y = data(:, 2);
z = data(:, 3);
u = data(:, 4);
v = data(:, 5);
w = data(:, 6);

% Create a grid for the data
[X, Y, Z] = meshgrid(unique(x), unique(y), unique(z));
U = griddata(x, y, z, u, X, Y, Z);
V = griddata(x, y, z, v, X, Y, Z);
W = griddata(x, y, z, w, X, Y, Z);

% Compute the magnitude of the velocity
speed = sqrt(U.^2 + V.^2 + W.^2);

% Avoid zero or very small values by setting a lower bound
speed(speed < 1e-6) = 1e-6;

%% 3D Visualization of the Flow
figure;
hold on;

% Slice views at specific planes
slice(X, Y, Z, speed, [], 0.5, []); % YZ plane (x = 0.5)
slice(X, Y, Z, speed, 0.5, [], []); % XZ plane (y = 0.5)
slice(X, Y, Z, speed, [], [], 0.5); % XY plane (z = 0.5)

% Add streamline for 3D visualization
seedX = linspace(min(x), max(x), 5);
seedY = linspace(min(y), max(y), 5);
seedZ = linspace(min(z), max(z), 5);
streamline(X, Y, Z, U, V, W, seedX, seedY, seedZ);

% Contours on the slices
contourslice(X, Y, Z, speed, [], 0.5, []); % Contour on YZ plane
contourslice(X, Y, Z, speed, 0.5, [], []); % Contour on XZ plane
contourslice(X, Y, Z, speed, [], [], 0.5); % Contour on XY plane

% Customize the visualization
shading interp; % Smooth color transitions
colormap(jet); % Color map for speed
colorbar; % Add color bar
title('3D Flow Field Visualization');
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
axis tight;
grid on;
view(3); % 3D view
hold off;

%% Plot XY Plane (z = 0.5)
figure;
% Extract the slice data for XY plane
[~, idx_z] = min(abs(Z(1, 1, :) - 0.5)); % Find the closest z-index to 0.5
speed_XY = squeeze(speed(:, :, idx_z)); % Convert to 2D matrix
U_XY = squeeze(U(:, :, idx_z)); % X-direction velocity (2D)
V_XY = squeeze(V(:, :, idx_z)); % Y-direction velocity (2D)

% Plot heatmap
imagesc(unique(x), unique(y), log10(speed_XY'));
set(gca, 'YDir', 'normal'); % Correct y-axis direction
hold on;

% Add streamlines
streamslice(unique(x), unique(y), U_XY', V_XY', 2);
title('XY Plane (z = 0.5)');
xlabel('X-axis');
ylabel('Y-axis');
colorbar;
colormap(turbo);
shading interp;
axis equal;
grid on;
hold off;