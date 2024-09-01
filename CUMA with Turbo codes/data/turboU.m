U_values = [10, 20, 30, 40, 50]; % Different values of U


% 绘制BER随L变化的图
figure;

% QPSK
plot(U_values, BER_results(1,:), '-s', 'LineWidth', 1.5, 'DisplayName', 'QPSK with Turbo');
hold on;
plot(U_values, BER_results(2,:), '--s', 'LineWidth', 1.5, 'DisplayName', 'QPSK without Turbo');

% 16QAM
plot(U_values, BER_results_16QAM(1,:), '-d', 'LineWidth', 1.5,  'DisplayName', '16QAM with Turbo');
plot(U_values, BER_results_16QAM(2,:), '--d', 'LineWidth', 1.5,  'DisplayName', '16QAM without Turbo');

% 64QAM
plot(U_values, BER_results_64QAM(1,:), '-^', 'LineWidth', 1.5, 'DisplayName', '64QAM with Turbo');
plot(U_values, BER_results_64QAM(2,:), '--^', 'LineWidth', 1.5, 'DisplayName', '64QAM without Turbo');

% 添加图例
legend;

% 添加标题和标签

xlabel('Number of Users (U)');
ylabel('Average BER');

% 显示网格
grid on;

hold off;
