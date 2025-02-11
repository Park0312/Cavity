% MATLAB Script to Automatically Load CSV Files and Save Plots as JPG
clc; clear; close all;

% 경로 설정
data_folder = './Animation_Result/';  % CSV 파일이 저장된 폴더
output_folder = './Plots/';           % JPG 파일을 저장할 폴더

% 폴더 생성 (존재하지 않으면 생성)
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% CSV 파일 목록 불러오기
csv_files = dir(fullfile(data_folder, '2D_results_9*.csv'));

% 각 CSV 파일에 대해 처리
for k = 1:length(csv_files)
    % 파일 경로 설정
    csv_file = fullfile(data_folder, csv_files(k).name);
    
    % 시간 값 추출 (파일 이름에서)
    time_str = extractBetween(csv_files(k).name, "2D_results_", "s.csv");
    time_value = str2double(time_str{1});

    % 데이터 불러오기
    data = readmatrix(csv_file);
    
    % 데이터 분해
    x = data(:, 1); % x-coordinates
    y = data(:, 2); % y-coordinates
    u = data(:, 3); % u-velocity
    v = data(:, 4); % v-velocity

    % 격자 크기 설정 (N+1 x N+1 형태)
    N = sqrt(length(x)) - 1;
    x = reshape(x, N+1, N+1);
    y = reshape(y, N+1, N+1);
    u = reshape(u, N+1, N+1);
    v = reshape(v, N+1, N+1);

    % 속도 크기 계산
    speed = sqrt(u.^2 + v.^2);

    % 📌 **속도 컨투어 플롯 & 벡터 그래프**
    figure('Visible', 'off'); % GUI 창을 띄우지 않고 처리
    contourf(x, y, speed, 20, 'LineStyle', 'none'); % 속도 컨투어 플롯
    colorbar;
    colormap("jet");
    hold on;
    quiver(x, y, u, v, 'k'); % 속도 벡터 (검은색 화살표)
    title(sprintf('Time = %.3fs: Speed Contour & Velocity Vectors', time_value));
    xlabel('X');
    ylabel('Y');
    axis equal;
    grid on;
    hold off;

    % 파일 저장 (JPG)
    saveas(gcf, fullfile(output_folder, sprintf('Speed_%.3fs.jpg', time_value)));

    % 📌 **y = 0.5에서 v-velocity 플롯**
    figure('Visible', 'off');
    y_target = 0.5;
    tolerance = 1e-6;
    indices_y = find(abs(y - y_target) < tolerance);
    x_at_y_half = x(indices_y);
    v_at_y_half = v(indices_y);
    plot(x_at_y_half, v_at_y_half, 'b-', 'LineWidth', 2);
    title(sprintf('Re=1000, V-velocity at y=H/2 (Time=%.3fs)', time_value));
    xlabel('x/L');
    ylabel('V/V0');
    grid on;
    saveas(gcf, fullfile(output_folder, sprintf('V_Velocity_%.3fs.jpg', time_value)));

    % 📌 **x = 0.5에서 u-velocity 플롯**
    figure('Visible', 'off');
    x_target = 0.5;
    indices_x = find(abs(x - x_target) < tolerance);
    y_at_x_half = y(indices_x);
    u_at_x_half = u(indices_x);
    plot(u_at_x_half, y_at_x_half, 'r-', 'LineWidth', 2);
    title(sprintf('Re=1000, U-velocity at x=H/2 (Time=%.3fs)', time_value));
    xlabel('U/U0');
    ylabel('y/L');
    grid on;
    saveas(gcf, fullfile(output_folder, sprintf('U_Velocity_%.3fs.jpg', time_value)));

    % 진행 상황 출력
    fprintf('✅ Processed and saved plots for: %s\n', csv_files(k).name);
end

fprintf('🎉 All plots have been saved in: %s\n', output_folder);