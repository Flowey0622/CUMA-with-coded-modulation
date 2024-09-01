U = 20; % the number of users
Nt = U; % number of transmit antennas (base station)
K = 0; % Rician factor for the channel model
L = 50; % Number of NLOS paths
lamda = 5; % Wavelength of radiation

% SNR values
snr_values = 0:5:30; % Range of SNR values in dB
numExperiments = 1000; % Number of experiments

% Port settings
N1 = 20; % Fixed port count
N2 = 20;
N = N1 * N2; % Total number of ports
modulation_schemes = {'BPSK', 'QPSK', '16-QAM', '64-QAM', '256-QAM'}; % Different modulation schemes
BER_results = zeros(length(modulation_schemes), length(snr_values)); % Store BER results

for m = 1:length(modulation_schemes)
    modulation_scheme = modulation_schemes{m};
    fprintf('Running simulation for modulation: %s\n', modulation_scheme);

    for s = 1:length(snr_values)
        snr = snr_values(s);
        
        % Initialize variables
        n1 = transpose((0:N1-1)./(N1-1).*10*lamda); % position vector n1
        n2 = transpose((0:N2-1)./(N2-1).*10*lamda); % position vector n2
        rho = 0.1; % Parameter to control the minimum acceptable level(0~1)

        K1 = cell(1, U);
        K2 = cell(1, U);
        H = cell(1, U);
        original_bits = cell(1, Nt);
        s_u = cell(1, Nt);
        demodulated_bits = cell(1, U);
        BER_users = zeros(U, numExperiments);

        for experiment = 1:numExperiments
            % Generate H matrix for each user
            for u = 1:U
                delta = 2 * pi * rand;
                theta = -pi + 2 * pi * rand;
                phi = -pi/2 + pi * rand;

                a1_rx = exp(-1i * 2 * pi .* n2 * sin(theta) * cos(phi));
                a2_rx = exp(-1i * 2 * pi .* n1 * sin(phi));
                a_rx = kron(a1_rx, a2_rx);
                n_tx = (0:Nt-1).';
                a_tx = exp(-1i * 2 * pi .* n_tx * sin(theta) * cos(phi));

                H_los = sqrt(K/(K+1)) * exp(-1i * delta) * (a_rx * a_tx');
                H_nlos = 0;

                for i = 1:L
                    kappa = sqrt(1/2) * complex(randn, randn);
                    theta = 2 * pi * rand;
                    phi = asin(2 * rand - 1);
                    a1_rx = exp(-1i * 2 * pi .* n2 * sin(theta) * cos(phi));
                    a2_rx = exp(-1i * 2 * pi .* n1 * sin(phi));
                    a_rx = kron(a1_rx, a2_rx);
                    a_tx = exp(-1i * 2 * pi .* n_tx * sin(theta) * cos(phi));
                    H_nlos = H_nlos + kappa * (a_rx * a_tx');
                end

                H_nlos = sqrt(1/(L*(K+1))) * H_nlos;
                H{u} = H_los + H_nlos;
            end

            % Generate a random orthonormal basis for the Nt x Nt complex space
            [Q, ~] = qr(complex(randn(Nt, Nt), randn(Nt, Nt)));
            for i = 1:Nt
                for j = 1:i-1
                    Q(i, :) = Q(i, :) - (Q(i, :) * Q(j, :)') * conj(Q(j, :));
                end
                Q(i, :) = Q(i, :) / norm(Q(i, :));
            end

            % Generate precoding vectors
            precoding_vectors = cell(1, U);
            g = cell(1, U);
            for u = 1:U
                precoding_vectors{u} = Q(:, u);
                g{u} = H{u} * precoding_vectors{u};
            end

            % Port selection based on real and imaginary parts
            for u = 1:U
                g_u = g{u};
                real_parts = real(g_u);
                imag_parts = imag(g_u);

                % Real part selection
                positive_indices1 = find(real_parts > 0);
                negative_indices1 = find(real_parts <= 0);

                if ~isempty(positive_indices1)
                    max_real_K1_plus = max(real_parts(positive_indices1));
                else
                    max_real_K1_plus = 0;
                end

                if ~isempty(negative_indices1)
                    min_real_K1_minus = min(real_parts(negative_indices1));
                else
                    min_real_K1_minus = 0;
                end

                K1_plus_selected = positive_indices1(real_parts(positive_indices1) >= rho * max_real_K1_plus);
                K1_minus_selected = negative_indices1(real_parts(negative_indices1) <= rho * min_real_K1_minus);

                sum_real_K1_plus = sum(real_parts(K1_plus_selected));
                sum_real_K1_minus = sum(real_parts(K1_minus_selected));

                if abs(sum_real_K1_plus) >= abs(sum_real_K1_minus)
                    K1{u} = K1_plus_selected;
                else
                    K1{u} = K1_minus_selected;
                end

                % Imaginary part selection
                positive_indices2 = find(imag_parts > 0);
                negative_indices2 = find(imag_parts <= 0);

                if ~isempty(positive_indices2)
                    max_imag_K2_plus = max(imag_parts(positive_indices2));
                else
                    max_imag_K2_plus = 0;
                end

                if ~isempty(negative_indices2)
                    min_imag_K2_minus = min(imag_parts(negative_indices2));
                else
                    min_imag_K2_minus = 0;
                end

                K2_plus_selected = positive_indices2(imag_parts(positive_indices2) >= rho * max_imag_K2_plus);
                K2_minus_selected = negative_indices2(imag_parts(negative_indices2) <= rho * min_imag_K2_minus);

                sum_imag_K2_plus = sum(imag_parts(K2_plus_selected));
                sum_imag_K2_minus = sum(imag_parts(K2_minus_selected));

                if abs(sum_imag_K2_plus) >= abs(sum_imag_K2_minus)
                    K2{u} = K2_plus_selected;
                else
                    K2{u} = K2_minus_selected;
                end
            end

            % Generate random binary data for each user and modulate
            for n = 1:Nt
                numBits = 0;
                numSymbols = 0;
                
                if strcmp(modulation_scheme, 'BPSK')
                    numBits = 100; % for BPSK, 1 bit per symbol
                    numSymbols = numBits;
                    original_bits{n} = randi([0 1], 1, numBits);
                    for k = 1:numSymbols
                        bit = original_bits{n}(k);
                        if bit == 0
                            s_u{n}(k) = 1;
                        else
                            s_u{n}(k) = -1;
                        end
                    end
                elseif strcmp(modulation_scheme, 'QPSK')
                    numBits = 200; % for QPSK, 2 bits per symbol
                    numSymbols = numBits / 2;
                    original_bits{n} = randi([0 1], 1, numBits);
                    s_u{n} = zeros(1, numSymbols);
                    for k = 1:numSymbols
                        bits_pair = original_bits{n}(2*k-1:2*k);
                        if isequal(bits_pair, [0 0])
                            s_u{n}(k) = 1 + 1i;
                        elseif isequal(bits_pair, [0 1])
                            s_u{n}(k) = 1 - 1i;
                        elseif isequal(bits_pair, [1 0])
                            s_u{n}(k) = -1 + 1i;
                        else
                            s_u{n}(k) = -1 - 1i;
                        end
                    end
                elseif strcmp(modulation_scheme, '16-QAM')
                    numBits = 400; % for 16-QAM, 4 bits per symbol
                    numSymbols = numBits / 4;
                    original_bits{n} = randi([0 1], 1, numBits);
                    s_u{n} = zeros(1, numSymbols);
                    for k = 1:numSymbols
                        bits = original_bits{n}(4*k-3:4*k);
                        if isequal(bits, [0 0 0 0])
                            s_u{n}(k) = 3 + 3i;
                        elseif isequal(bits, [0 0 0 1])
                            s_u{n}(k) = 3 + 1i;
                        elseif isequal(bits, [0 0 1 0])
                            s_u{n}(k) = 3 - 3i;
                        elseif isequal(bits, [0 0 1 1])
                            s_u{n}(k) = 3 - 1i;
                        elseif isequal(bits, [0 1 0 0])
                            s_u{n}(k) = 1 + 3i;
                        elseif isequal(bits, [0 1 0 1])
                            s_u{n}(k) = 1 + 1i;
                        elseif isequal(bits, [0 1 1 0])
                            s_u{n}(k) = 1 - 3i;
                        elseif isequal(bits, [0 1 1 1])
                            s_u{n}(k) = 1 - 1i;
                        elseif isequal(bits, [1 0 0 0])
                            s_u{n}(k) = -1 + 3i;
                        elseif isequal(bits, [1 0 0 1])
                            s_u{n}(k) = -1 + 1i;
                        elseif isequal(bits, [1 0 1 0])
                            s_u{n}(k) = -1 - 3i;
                        elseif isequal(bits, [1 0 1 1])
                            s_u{n}(k) = -1 - 1i;
                        elseif isequal(bits, [1 1 0 0])
                            s_u{n}(k) = -3 + 3i;
                        elseif isequal(bits, [1 1 0 1])
                            s_u{n}(k) = -3 + 1i;
                        elseif isequal(bits, [1 1 1 0])
                            s_u{n}(k) = -3 - 3i;
                        else
                            s_u{n}(k) = -3 - 1i;
                        end
                    end
                elseif strcmp(modulation_scheme, '64-QAM')
                    numBits = 600; % for 64-QAM, 6 bits per symbol
                    numSymbols = numBits / 6;
                    original_bits{n} = randi([0 1], 1, numBits);
                    s_u{n} = zeros(1, numSymbols);
                    for k = 1:numSymbols
                        bits = original_bits{n}(6*k-5:6*k);
                        real_part = 0;
                        imag_part = 0;
                        if isequal(bits(1:3), [0 0 0])
                            real_part = 7;
                        elseif isequal(bits(1:3), [0 0 1])
                            real_part = 5;
                        elseif isequal(bits(1:3), [0 1 0])
                            real_part = 3;
                        elseif isequal(bits(1:3), [0 1 1])
                            real_part = 1;
                        elseif isequal(bits(1:3), [1 1 0])
                            real_part = -3;
                        elseif isequal(bits(1:3), [1 1 1])
                            real_part = -1;
                        elseif isequal(bits(1:3), [1 0 0])
                            real_part = -7;
                        else
                            real_part = -5;
                        end
                        if isequal(bits(4:6), [0 0 0])
                            imag_part = 7;
                        elseif isequal(bits(4:6), [0 0 1])
                            imag_part = 5;
                        elseif isequal(bits(4:6), [0 1 0])
                            imag_part = 3;
                        elseif isequal(bits(4:6), [0 1 1])
                            imag_part = 1;
                        elseif isequal(bits(4:6), [1 1 0])
                            imag_part = -3;
                        elseif isequal(bits(4:6), [1 1 1])
                            imag_part = -1;
                        elseif isequal(bits(4:6), [1 0 0])
                            imag_part = -7;
                        else
                            imag_part = -5;
                        end
                        s_u{n}(k) = real_part + 1i * imag_part;
                    end
                elseif strcmp(modulation_scheme, '256-QAM')
                    numBits = 800; % for 256-QAM, 8 bits per symbol
                    numSymbols = numBits / 8;
                    original_bits{n} = randi([0 1], 1, numBits);
                    s_u{n} = zeros(1, numSymbols);
                    for k = 1:numSymbols
                        bits = original_bits{n}(8*k-7:8*k);
                        real_part = 0;
                        imag_part = 0;
                        if isequal(bits(1:4), [0 0 0 0])
                            real_part = 15;
                        elseif isequal(bits(1:4), [0 0 0 1])
                            real_part = 13;
                        elseif isequal(bits(1:4), [0 0 1 0])
                            real_part = 11;
                        elseif isequal(bits(1:4), [0 0 1 1])
                            real_part = 9;
                        elseif isequal(bits(1:4), [0 1 0 0])
                            real_part = 7;
                        elseif isequal(bits(1:4), [0 1 0 1])
                            real_part = 5;
                        elseif isequal(bits(1:4), [0 1 1 0])
                            real_part = 3;
                        elseif isequal(bits(1:4), [0 1 1 1])
                            real_part = 1;
                        elseif isequal(bits(1:4), [1 0 0 0])
                            real_part = -15;
                        elseif isequal(bits(1:4), [1 0 0 1])
                            real_part = -13;
                        elseif isequal(bits(1:4), [1 0 1 0])
                            real_part = -11;
                        elseif isequal(bits(1:4), [1 0 1 1])
                            real_part = -9;
                        elseif isequal(bits(1:4), [1 1 0 0])
                            real_part = -7;
                        elseif isequal(bits(1:4), [1 1 0 1])
                            real_part = -5;
                        elseif isequal(bits(1:4), [1 1 1 0])
                            real_part = -3;
                        else
                            real_part = -1;
                        end
                        if isequal(bits(5:8), [0 0 0 0])
                            imag_part = 15;
                        elseif isequal(bits(5:8), [0 0 0 1])
                            imag_part = 13;
                        elseif isequal(bits(5:8), [0 0 1 0])
                            imag_part = 11;
                        elseif isequal(bits(5:8), [0 0 1 1])
                            imag_part = 9;
                        elseif isequal(bits(5:8), [0 1 0 0])
                            imag_part = 7;
                        elseif isequal(bits(5:8), [0 1 0 1])
                            imag_part = 5;
                        elseif isequal(bits(5:8), [0 1 1 0])
                            imag_part = 3;
                        elseif isequal(bits(5:8), [0 1 1 1])
                            imag_part = 1;
                        elseif isequal(bits(5:8), [1 0 0 0])
                            imag_part = -15;
                        elseif isequal(bits(5:8), [1 0 0 1])
                            imag_part = -13;
                        elseif isequal(bits(5:8), [1 0 1 0])
                            imag_part = -11;
                        elseif isequal(bits(5:8), [1 0 1 1])
                            imag_part = -9;
                        elseif isequal(bits(5:8), [1 1 0 0])
                            imag_part = -7;
                        elseif isequal(bits(5:8), [1 1 0 1])
                            imag_part = -5;
                        elseif isequal(bits(5:8), [1 1 1 0])
                            imag_part = -3;
                        else
                            imag_part = -1;
                        end
                        s_u{n}(k) = real_part + 1i * imag_part;
                    end
                end
            end

            eta = cell(1, U);
            sigma = 0.1; 
            signal_power = 1;
            snr_linear = 10^(snr / 10); 
            noise_power = signal_power / snr_linear; 

            for u = 1:U
                real_noise = sqrt(noise_power / 2) * randn(N, numSymbols); 
                imag_noise = sqrt(noise_power / 2) * randn(N, numSymbols); 
                eta{u} = real_noise + 1i * imag_noise; 
            end
            total_interference = cell(1, U);

            % Simulate transmission through the channel with AWGN
            for u = 1:U
                g_u = g{u};
                g_u_u1 = g_u(K1{u});
                g_u_u2 = g_u(K2{u});

                g_u_u1_real = real(g_u_u1);
                g_u_u1_imag = imag(g_u_u1);
                g_u_u2_real = real(g_u_u2);
                g_u_u2_imag = imag(g_u_u2);

                s_u_I = real(s_u{u});
                s_u_Q = imag(s_u{u});

                A = sum(g_u_u1_real);
                a = sum(g_u_u1_imag);
                b = sum(g_u_u2_real);
                B = sum(g_u_u2_imag);

                r_I1 = A * s_u_I - a * s_u_Q;
                r_Q1 = a * s_u_I + A * s_u_Q;
                r_I2 = b * s_u_I - B * s_u_Q;
                r_Q2 = B * s_u_I + b * s_u_Q;

                % Total interference from other users
                total_interference{u} = eta{u};

                for other_u = 1:U
                    if other_u ~= u
                        g_u_other = H{u} * precoding_vectors{other_u};
                        s_other = s_u{other_u};
                        interference = g_u_other * s_other; % Cumulative interference
                        total_interference{u} = total_interference{u} + interference;
                    end
                end

                noise_I1 = sum(real(total_interference{u}(K1{u})));
                noise_Q1 = sum(imag(total_interference{u}(K1{u})));
                noise_I2 = sum(real(total_interference{u}(K2{u})));
                noise_Q2 = sum(imag(total_interference{u}(K2{u})));

                received_signal_I1{u} = r_I1 + noise_I1;
                received_signal_Q1{u} = r_Q1 + noise_Q1;
                received_signal_I2{u} = r_I2 + noise_I2;
                received_signal_Q2{u} = r_Q2 + noise_Q2;

                r_u{u} = [received_signal_I1{u}; 
                          received_signal_Q1{u};   
                          received_signal_I2{u}; 
                          received_signal_Q2{u}];  % 4*numsymbols

                G_u_u = [A, -a;
                         a,  A;
                         b, -B;
                         B,  b]; % 4*2

                estimated_symbols{u} = G_u_u \ r_u{u};

                % Demodulate
                if strcmp(modulation_scheme, 'BPSK')
                    demodulated_bits{u} = zeros(1,  numSymbols);
                    for h =1: numSymbols
                        real_part = estimated_symbols{u}(1, h);
                        if real_part > 0
                            demodulated_bits{u}(h) = 0;
                        else
                            demodulated_bits{u}(h) = 1;
                        end
                    end
                    
                elseif strcmp(modulation_scheme, 'QPSK')
                    demodulated_bits{u} = zeros(1, 2 * numSymbols);
                    for h = 1:numSymbols
                        real_part = estimated_symbols{u}(1, h);
                        imag_part = imag(estimated_symbols{u}(2, h));
                        if real_part > 0 && imag_part > 0
                            demodulated_bits{u}(2*h-1:2*h) = [0 0]; % 1+1i
                        elseif real_part > 0 && imag_part < 0
                            demodulated_bits{u}(2*h-1:2*h) = [0 1]; % 1-1i
                        elseif real_part < 0 && imag_part > 0
                            demodulated_bits{u}(2*h-1:2*h) = [1 0]; % -1+1i
                        else
                            demodulated_bits{u}(2*h-1:2*h) = [1 1]; % -1-1i
                        end
                    end
                elseif strcmp(modulation_scheme, '16-QAM')
                    demodulated_bits{u} = zeros(1, 4 * numSymbols);
                    for h = 1:numSymbols
                        real_part = estimated_symbols{u}(1, h);
                        imag_part = estimated_symbols{u}(2, h);
                        if real_part > 2
                            real_bits = [0 0];
                        elseif real_part > 0
                            real_bits = [0 1];
                        elseif real_part > -2
                            real_bits = [1 0];
                        else
                            real_bits = [1 1];
                        end

                        if imag_part > 2
                            imag_bits = [0 0];
                        elseif imag_part > 0
                            imag_bits = [0 1];
                        elseif imag_part > -2
                            imag_bits = [1 0];
                        else
                            imag_bits = [1 1];
                        end

                        demodulated_bits{u}(4*h-3:4*h) = [real_bits, imag_bits];
                    end
                elseif strcmp(modulation_scheme, '64-QAM')
                    demodulated_bits{u} = zeros(1, 6 * numSymbols);
                    for h = 1:numSymbols
                        real_part = estimated_symbols{u}(1, h);
                        imag_part = estimated_symbols{u}(2, h);
                        if real_part > 6
                            real_bits = [0 0 0];
                        elseif real_part > 4
                            real_bits = [0 0 1];
                        elseif real_part > 2
                            real_bits = [0 1 0];
                        elseif real_part > 0
                            real_bits = [0 1 1];
                        elseif real_part > -2
                            real_bits = [1 1 1];
                        elseif real_part > -4
                            real_bits = [1 1 0];
                        elseif real_part > -6
                            real_bits = [1 0 1];
                        else
                            real_bits = [1 0 0];
                        end

                        if imag_part > 6
                            imag_bits = [0 0 0];
                        elseif imag_part > 4
                            imag_bits = [0 0 1];
                        elseif imag_part > 2
                            imag_bits = [0 1 0];
                        elseif imag_part > 0
                            imag_bits = [0 1 1];
                        elseif imag_part > -2
                            imag_bits = [1 1 1];
                        elseif imag_part > -4
                            imag_bits = [1 1 0];
                        elseif imag_part > -6
                            imag_bits = [1 0 1];
                        else
                            imag_bits = [1 0 0];
                        end

                        demodulated_bits{u}(6*h-5:6*h) = [real_bits, imag_bits];
                    end
                elseif strcmp(modulation_scheme, '256-QAM')
                    demodulated_bits{u} = zeros(1, 8 * numSymbols);
                    for h = 1:numSymbols
                        real_part = estimated_symbols{u}(1, h);
                        imag_part = estimated_symbols{u}(2, h);
                        if real_part > 14
                            real_bits = [0 0 0 0];
                        elseif real_part > 12
                            real_bits = [0 0 0 1];
                        elseif real_part > 10
                            real_bits = [0 0 1 0];
                        elseif real_part > 8
                            real_bits = [0 0 1 1];
                        elseif real_part > 6
                            real_bits = [0 1 0 0];
                        elseif real_part > 4
                            real_bits = [0 1 0 1];
                        elseif real_part > 2
                            real_bits = [0 1 1 0];
                        elseif real_part > 0
                            real_bits = [0 1 1 1];
                        elseif real_part > -2
                            real_bits = [1 1 1 1];
                        elseif real_part > -4
                            real_bits = [1 1 1 0];
                        elseif real_part > -6
                            real_bits = [1 1 0 1];
                        elseif real_part > -8
                            real_bits = [1 1 0 0];
                        elseif real_part > -10
                            real_bits = [1 0 1 1];
                        elseif real_part > -12
                            real_bits = [1 0 1 0];
                        elseif real_part > -14
                            real_bits = [1 0 0 1];
                        else
                            real_bits = [1 0 0 0];
                        end

                        if imag_part > 14
                            imag_bits = [0 0 0 0];
                        elseif imag_part > 12
                            imag_bits = [0 0 0 1];
                        elseif imag_part > 10
                            imag_bits = [0 0 1 0];
                        elseif imag_part > 8
                            imag_bits = [0 0 1 1];
                        elseif imag_part > 6
                            imag_bits = [0 1 0 0];
                        elseif imag_part > 4
                            imag_bits = [0 1 0 1];
                        elseif imag_part > 2
                            imag_bits = [0 1 1 0];
                        elseif imag_part > 0
                            imag_bits = [0 1 1 1];
                        elseif imag_part > -2
                            imag_bits = [1 1 1 1];
                        elseif imag_part > -4
                            imag_bits = [1 1 1 0];
                        elseif imag_part > -6
                            imag_bits = [1 1 0 1];
                        elseif imag_part > -8
                            imag_bits = [1 1 0 0];
                        elseif imag_part > -10
                            imag_bits = [1 0 1 1];
                        elseif imag_part > -12
                            imag_bits = [1 0 1 0];
                        elseif imag_part > -14
                            imag_bits = [1 0 0 1];
                        else
                            imag_bits = [1 0 0 0];
                        end

                        demodulated_bits{u}(8*h-7:8*h) = [real_bits, imag_bits];
                    end
                end

                % Calculate number of bit errors
                num_errors = sum(original_bits{u} ~= demodulated_bits{u});
                % Calculate BER
                BER_users(u, experiment) = num_errors / numBits;
            end
        end

        % Calculate average BER for this SNR and modulation scheme
        average_BER_users = mean(BER_users, 2);
        average_BER = mean(average_BER_users);
        BER_results(m, s) = average_BER;
    end
end

% Plotting the results
figure;
hold on;
markers = {'-o', '-s', '-d', '-^', '-v'};
for m = 1:length(modulation_schemes)
    plot(snr_values, BER_results(m, :), markers{m}, 'LineWidth', 2, 'DisplayName', modulation_schemes{m});
end
set(gca, 'YScale', 'log');
xlabel('SNR (dB)');
ylabel('Average BER');
legend('show');
grid on;
hold off;
