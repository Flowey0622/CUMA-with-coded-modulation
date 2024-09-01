% 定义NLOS路径数量
L_values = [10, 30, 50, 70, 100];


% 绘制BER随L变化的图
figure;

% QPSK
plot(L_values, BER_results_QPSK(1,:), '-s', 'LineWidth', 1.5, 'DisplayName', 'QPSK with LDPC');
hold on;
plot(L_values, BER_results_QPSK(2,:), '--s', 'LineWidth', 1.5, 'DisplayName', 'QPSK without LDPC');

% 16QAM
plot(L_values, BER_results_16QAM(1,:), '-d', 'LineWidth', 1.5,  'DisplayName', '16QAM with LDPC');
plot(L_values, BER_results_16QAM(2,:), '--d', 'LineWidth', 1.5,  'DisplayName', '16QAM without LDPC');

% 64QAM
plot(L_values, BER_results_64QAM(1,:), '-^', 'LineWidth', 1.5, 'DisplayName', '64QAM with LDPC');
plot(L_values, BER_results_64QAM(2,:), '--^', 'LineWidth', 1.5, 'DisplayName', '64QAM without LDPC');

% 添加图例
legend;

% 添加标题和标签

xlabel('Number of NLOS Paths (L)');
ylabel('Average BER');

% 显示网格
grid on;

hold off;
