port_counts = [10, 20, 30, 40, 50];

figure;

% QPSK
plot(port_counts.^2, BER_results_QPSK(1,:), '-s', 'LineWidth', 1.5, 'DisplayName', 'QPSK with LDPC');
hold on;
plot(port_counts.^2, BER_results_QPSK(2,:), '--s', 'LineWidth', 1.5, 'DisplayName', 'QPSK without LDPC');

% 16QAM
plot(port_counts.^2, BER_results_16QAM(1,:), '-d', 'LineWidth', 1.5,  'DisplayName', '16QAM with LDPC');
plot(port_counts.^2, BER_results_16QAM(2,:), '--d', 'LineWidth', 1.5,  'DisplayName', '16QAM without LDPC');

% 64QAM
plot(port_counts.^2, BER_results_64QAM(1,:), '-^', 'LineWidth', 1.5, 'DisplayName', '64QAM with LDPC');
plot(port_counts.^2, BER_results_64QAM(2,:), '--^', 'LineWidth', 1.5, 'DisplayName', '64QAM without LDPC');


legend;
xlabel('Number of Ports (N)');
ylabel('Average BER');
grid on;
hold off;
