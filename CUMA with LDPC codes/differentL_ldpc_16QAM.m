P = [16 17 22 24  9  3 14 -1  4  2  7 -1 26 -1  2 -1 21 -1  1  0 -1 -1 -1 -1
     25 12 12  3  3 26  6 21 -1 15 22 -1 15 -1  4 -1 -1 16 -1  0  0 -1 -1 -1
     25 18 26 16 22 23  9 -1  0 -1  4 -1  4 -1  8 23 11 -1 -1 -1  0  0 -1 -1
      9  7  0  1 17 -1 -1  7  3 -1  3 23 -1 16 -1 -1 21 -1  0 -1 -1  0  0 -1
     24  5 26  7  1 -1 -1 15 24 15 -1  8 -1 13 -1 13 -1 11 -1 -1 -1 -1  0  0
      2  2 19 14 24  1 15 19 -1 21 -1  2 -1 24 -1  3 -1  2  1 -1 -1 -1 -1  0];
blockSize = 27;
pcmatrix = ldpcQuasiCyclicMatrix(blockSize,P);

encodercfg = ldpcEncoderConfig(pcmatrix);
decodercfg = ldpcDecoderConfig(pcmatrix);
maxnumiter = 100;
M=16;

U=20; % the number of user
Nt = U; % number of transmit antennas (base station)

%finite scattering channel model
%environment parameters
K = 0;
% NLOS path settings
L_values = [10, 30, 50, 70, 100]; % Different NLOS paths
BER_results_16QAM = zeros(2, length(L_values)); % Store BER results

% fluid antenna parameters: only Rx parameters are FAS
lamda = 5; % wavelength of radiation
%physical dimensions of the space
W1 = 10 *lamda; % 3 lambda
W2 = 10 *lamda; % 1.6 lambda
W = W1*W2; %2D FAS for each UE

% number of ports in each dimension
N1=20;
N2=20;%the number of ports at each user
N=N1*N2; % total number 
n1 = transpose((0:N1-1)./(N1-1).*W1); % position vector
n2 = transpose((0:N2-1)./(N2-1).*W2);

% Select K ports that satify the condition
% seleted for activation
rho = 0.01; % Parameter to control the minimum acceptable level(0~1)

% Simulation parameters
% Set the message length based on the encoder configuration
messageLength = encodercfg.NumInformationBits;
numBits = messageLength; 
snr = 30; % Signal-to-noise ratio in dB

%available ports to obtain the received signal
K1 = cell(1, U); % more focus on in-phase component of the channel
K2 = cell(1, U);  % more focus on quadrature component of the channel
% Initialize the port sets
K1_plus = cell(1,Nt);  % K+
K1_minus = cell(1,Nt);  % K-
K2_plus = cell(1,Nt);  % K+
K2_minus = cell(1,Nt);  % K-
sum_real_K1_plus = zeros(Nt, 1);
sum_real_K1_minus = zeros(Nt, 1);
sum_imag_K2_plus = zeros(Nt, 1);
sum_imag_K2_minus = zeros(Nt, 1);
% Initialize channel matrices for each user
H = cell(1, U);
% different infoSymbols for each user
infoSymbols = cell(1, Nt);
s_u_ldpc = cell(1,Nt);
s_u_ldpc2 = cell(1,Nt);
s_u = cell(1,Nt);
original_bits = cell(1,Nt);
encoded_bits = cell(1,Nt);
coded_bits = cell(1,Nt);
s_u_turbo = cell(1,Nt);
padded_bits = cell(1,Nt);
% Initialize storage for QPSK signals and received data
r_u = cell(1, U);
r_u_ldpc = cell(1, U);
r_u_turbo = cell(1, U);
R_U =cell(1,U);
received_signal_I1_ldpc = cell(1, U);
received_signal_Q1_ldpc = cell(1, U);
received_signal_I2_ldpc = cell(1, U);
received_signal_Q2_ldpc = cell(1, U);
received_signal_I1_turbo = cell(1, U);
received_signal_Q1_turbo = cell(1, U);
received_signal_I2_turbo = cell(1, U);
received_signal_Q2_turbo = cell(1, U);
received_signal_I1 = cell(1, U);
received_signal_Q1 = cell(1, U);
received_signal_I2 = cell(1, U);
received_signal_Q2 = cell(1, U);

demodulated_bits = cell(1, U);
decoded_bits_ldpc= cell(1, U);
decoded_bits_turbo= cell(1, U);
% decode symbols
estimated_symbols = cell(1, U);
estimated_symbols_ldpc = cell(1, U);
estimated_symbols_turbo = cell(1, U);
decoded_signal = cell(1, U);
BER_ldpc = zeros(1, U);
BER_turbo = zeros(1, U);
BER = zeros(1, U);
numExperiments = 1000; % Number of experiments
BER_user1 = zeros(1, numExperiments); % Array to store BER for user 1
BER_users_ldpc = zeros(U, numExperiments);
BER_users_turbo = zeros(U, numExperiments);
BER_users = zeros(U, numExperiments);


for l = 1:length(L_values)
    L = L_values(l);
    for experiment = 1:numExperiments
        for u = 1:U
            %los, azimuth, elevation
            delta = 2*pi*rand;
            theta = -pi+2*pi*rand;
            phi = -pi/2+pi*rand;
            %steering vector for the receiver
            a1_rx = [exp(-1i*2*pi.*n2*sin(theta)*cos(phi))];
            a2_rx = [exp(-1i*2*pi.*n1*sin(phi))];
            a_rx = kron(a1_rx,a2_rx);
            %steering vector for the transmitter
            n_tx = (0:Nt-1).'; % assuming linear array for Tx antennas
            a_tx = exp(-1i*2*pi.*n_tx*sin(theta)*cos(phi));
            
            %channel
            H_los = sqrt(K/(K+1))*exp(-1i*delta)*(a_rx*a_tx');
            H_nlos = 0;
            for i =1:L
                kappa = sqrt(1/2)*complex(randn,randn); %nlos coefficient
                theta = 2*pi*rand;
                phi = asin(2*rand-1); %azimuth,elevation
                a1_rx = [exp(-1i*2*pi.*n2*sin(theta)*cos(phi))];
                a2_rx = [exp(-1i*2*pi.*n1*sin(phi))];
                a_rx = kron(a1_rx,a2_rx);
                a_tx = exp(-1i*2*pi.*n_tx*sin(theta)*cos(phi));
                
                H_nlos = H_nlos +kappa*(a_rx * a_tx');
            end
            H_nlos = sqrt(1/(L*(K+1)))*H_nlos;
            H{u} = (H_los+H_nlos);
        end
        
        %for u = 1:U
         %   fprintf('Channel State Information (CSI) matrix H for user %d:\n', u);
          %  disp(H{u});
        %end
        
        % Generate a random orthonormal basis for the Nt x Nt complex space
        [Q, ~] = qr(complex(randn(Nt, Nt), randn(Nt, Nt))); 

        for i = 1:Nt
            for j = 1:i-1
                Q(i, :) = Q(i, :) - (Q(i, :) * Q(j, :)') * conj(Q(j, :));
            end
            Q(i, :) = Q(i, :) / norm(Q(i, :));
        end
        
        orthogonal_check = Q * Q';
        identity_matrix = eye(Nt);
        
        if norm(orthogonal_check - identity_matrix) < 1e-10
            disp('Q is an orthogonal matrix.');
        else
            disp('Q is not an orthogonal matrix.');
        end
        
        for i = 1:Nt
            for j = 1:Nt
                if i ~= j
                    dot_product = Q(i, :) * Q(j, :)';
                    %fprintf('Dot product of row %d and row %d: %f + %fi\n', i, j, real(dot_product), imag(dot_product));
                end
            end
        end
        
        % Generate random orthonormal precoding vector for each user
        precoding_vectors = cell(1, U);
        g = cell(1, U);
        for u = 1:U
            precoding_vectors{u} = 0.2*Q(:, u);
            
            % Calculate g for the current user
            g{u} = H{u} * precoding_vectors{u};
        end
        
        %K1 selection focusing on real parts
        %K2 selection focusing on imag parts
        for u = 1:U
            % Get the real and imag part of the channel matrix for the user
            g_u = g{u};
            real_parts = real(g_u);
            imag_parts = imag(g_u);
        
            % Initialize the selected ports list
            K1{u} = [];
            K2{u} = [];
         
            % Split ports into K_plus and K_minus based on the sign of their real parts
            positive_indices1 = find(real_parts > 0);
            negative_indices1 = find(real_parts <= 0);
        
            K1_plus{u} = positive_indices1;
            K1_minus{u} = negative_indices1;
            
            % To reduce the value of K
            K1_plus_selected = [];
            K1_minus_selected = [];
        
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
        
            % Select ports based on the given condition for K1
            K1_plus_selected = positive_indices1(real_parts(positive_indices1) >= rho * max_real_K1_plus);
            K1_minus_selected = negative_indices1(real_parts(negative_indices1) <= rho * min_real_K1_minus);
            
            % Calculate the sum of real parts for both sets for K1
            sum_real_K1_plus(u) = sum(real(g_u(K1_plus_selected)));
            sum_real_K1_minus(u) = sum(real(g_u(K1_minus_selected)));
            
            % Select the set with the larger absolute sum of real parts
            if abs(sum_real_K1_plus(u)) >= abs(sum_real_K1_minus(u))
                K1{u} = K1_plus_selected;
            else
                K1{u} = K1_minus_selected;
            end
        
            % Check that all selected ports have the same sign
            if ~isempty(K1{u})
                    selected_sign_K1 = sign(real_parts(K1{u}(1)));
                    for port = K1{u}
                        if sign(real_parts(port)) ~= selected_sign_K1
                            error('Selected ports do not have the same sign in real part.');
                        end
                    end
            end
        
            % Split ports into K_plus and K_minus based on the sign of their imag parts
            positive_indices2 = find(imag_parts > 0);
            negative_indices2 = find(imag_parts <= 0);
        
            K2_plus{u} = positive_indices2;
            K2_minus{u} = negative_indices2;
            
            % To reduce the value of K
            K2_plus_selected = [];
            K2_minus_selected = [];
        
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
        
            % Select ports based on the given condition for K2
            K2_plus_selected = positive_indices2(imag_parts(positive_indices2) >= rho * max_imag_K2_plus);
            K2_minus_selected = negative_indices2(imag_parts(negative_indices2) <= rho * min_imag_K2_minus);
            
            % Calculate the sum of imag parts for both sets for K2
            sum_imag_K2_plus(u) = sum(imag(g_u(K2_plus_selected)));
            sum_imag_K2_minus(u) = sum(imag(g_u(K2_minus_selected)));
            
            % Select the set with the larger absolute sum of imag parts
            if abs(sum_imag_K2_plus(u)) >= abs(sum_imag_K2_minus(u))
                K2{u} = K2_plus_selected;
            else
                K2{u} = K2_minus_selected;
            end
        
            % Check that all selected ports have the same sign
            if ~isempty(K2{u})
                    selected_sign_K2 = sign(imag_parts(K2{u}(1)));
                    for port = K2{u}
                        if sign(imag_parts(port)) ~= selected_sign_K2
                            error('Selected ports do not have the same sign in imag part.');
                        end
                    end
            end
        
            %fprintf('Selected K1 ports for user %d:\n', u);
            %disp(K1{u});
            %fprintf('Selected K2 ports for user %d:\n', u);
            %disp(K2{u});
        end
        
        
        % Generate random binary data and simulate downlink transmission for each user
        for n=1:Nt   
            original_bits{n} = randi([0 1], numBits,1); 
            
            % LDPC encode the data
            encoded_bits{n} = ldpcEncode(original_bits{n}, encodercfg);
    
            s_u_ldpc{n}= qammod(encoded_bits{n},M,InputType='bit');
            s_u{n} = qammod([original_bits{n}; zeros(2, 1)],M,InputType='bit');
            
            % ldpc
            numEncodedBits = length(encoded_bits{n});
            numEncodedSymbols = numEncodedBits / 4;
     
            numSymbols =length([original_bits{n}; zeros(2, 1)])/4;
           
        end
    
        eta_ldpc = cell(1, U);
        eta = cell(1, U);
        noise = cell(1,U);
        sigma = 0.1; 
        signal_power = 1; 
        snr_linear = 10^(snr / 10);
        noise_power = signal_power / snr_linear; 
        for u = 1:U
            real_noise = sqrt(noise_power / 2) * randn(N, numEncodedSymbols); 
            imag_noise = sqrt(noise_power / 2) * randn(N, numEncodedSymbols); 
            eta_ldpc{u} = real_noise + 1i * imag_noise;       
            eta{u} = eta_ldpc{u}(:, 1:numSymbols);
        end
        total_interference = cell(1,U);
        total_interference_ldpc = cell(1,U);
        
        % Simulate transmission through the channel with AWGN
        for u = 1:U
            g_u = g{u};
            % User's own signal contribution with K1 ports
            g_u_u1 = g_u(K1{u});
            g_u_u2 = g_u(K2{u});
        
            g_u_u1_real = real(g_u_u1);
            g_u_u1_imag = imag(g_u_u1);
            g_u_u2_real = real(g_u_u2);
            g_u_u2_imag = imag(g_u_u2);
        
            s_u_I_ldpc = real(s_u_ldpc{u}); % 1*100
            s_u_Q_ldpc = imag(s_u_ldpc{u});
            s_u_I = real(s_u{u}); % 1*100
            s_u_Q = imag(s_u{u});
        
            A = sum(g_u_u1_real);  
            a = sum(g_u_u1_imag);
            b = sum(g_u_u2_real);  
            B = sum(g_u_u2_imag);
        
            r_I1_ldpc = A * s_u_I_ldpc - a * s_u_Q_ldpc; % 100*1
            r_Q1_ldpc = a * s_u_I_ldpc + A * s_u_Q_ldpc;
            r_I2_ldpc = b * s_u_I_ldpc - B * s_u_Q_ldpc;
            r_Q2_ldpc = B * s_u_I_ldpc + b * s_u_Q_ldpc;
    
            r_I1 = A * s_u_I - a * s_u_Q; % 100*1
            r_Q1 = a * s_u_I + A * s_u_Q;
            r_I2 = b * s_u_I - B * s_u_Q;
            r_Q2 = B * s_u_I + b * s_u_Q;
        
            % Total interference from other users
            total_interference_ldpc{u}= eta_ldpc{u};
            total_interference{u}= eta{u};
        
            for other_u = 1:U
                if other_u ~= u
                    g_u_other = H{u} * precoding_vectors{other_u}; % 100*1
                    s_other_ldpc = s_u_ldpc{other_u};
                    interference = g_u_other * s_other_ldpc.'; % Cumulative interference
                    total_interference_ldpc{u} = total_interference_ldpc{u} + interference;
                end
            end
    
            for other_u = 1:U
                if other_u ~= u
                    g_u_other = H{u} * precoding_vectors{other_u}; % 100*1
                    s_other = s_u{other_u};
                    interference = g_u_other * s_other.'; % Cumulative interference
                    total_interference{u} = total_interference{u} +  interference;
                end
            end
        
            noise_I1_ldpc =sum(real(total_interference_ldpc{u}(K1{u},:)));
            noise_Q1_ldpc =sum(imag(total_interference_ldpc{u}(K1{u},:)));
            noise_I2_ldpc =sum(real(total_interference_ldpc{u}(K2{u},:)));
            noise_Q2_ldpc =sum(imag(total_interference_ldpc{u}(K2{u},:)));
    
            noise_I1 =sum(real(total_interference{u}(K1{u},:)));
            noise_Q1 =sum(imag(total_interference{u}(K1{u},:)));
            noise_I2 =sum(real(total_interference{u}(K2{u},:)));
            noise_Q2 =sum(imag(total_interference{u}(K2{u},:)));
    
            N0 = mean(mean(abs(total_interference_ldpc{u}).^2));

            received_signal_I1_ldpc{u} = r_I1_ldpc + noise_I1_ldpc.';
            received_signal_Q1_ldpc{u} = r_Q1_ldpc + noise_Q1_ldpc.';
            received_signal_I2_ldpc{u} = r_I2_ldpc + noise_I2_ldpc.';
            received_signal_Q2_ldpc{u} = r_Q2_ldpc + noise_Q2_ldpc.';
    
            received_signal_I1{u} = r_I1 + noise_I1.';
            received_signal_Q1{u} = r_Q1 + noise_Q1.';
            received_signal_I2{u} = r_I2 + noise_I2.';
            received_signal_Q2{u} = r_Q2 + noise_Q2.';
            
            r_u_ldpc{u} = [received_signal_I1_ldpc{u}.'; 
                           received_signal_Q1_ldpc{u}.';   
                           received_signal_I2_ldpc{u}.'; 
                           received_signal_Q2_ldpc{u}.'];  % 4*numsymbols
    
            r_u{u} = [received_signal_I1{u}.'; 
                      received_signal_Q1{u}.';   
                      received_signal_I2{u}.'; 
                      received_signal_Q2{u}.'];  % 4*numsymbols
         
            G_u_u = [A, -a;
                     a,  A;
                     b, -B;
                     B,  b]; % 4*2
            
            estimated_symbols{u}= G_u_u \ r_u{u};
            estimated_symbols_ldpc{u} = G_u_u \ r_u_ldpc{u};
          
            estimated_symbols{u}= estimated_symbols{u}(1, :)+ 1i *estimated_symbols{u}(2,:);
            estimated_symbols_ldpc{u}= estimated_symbols_ldpc{u}(1, :)+ 1i *estimated_symbols_ldpc{u}(2,:);
           
            demodulatedLLR_ldpc = qamdemod(estimated_symbols_ldpc{u}.', M, 'OutputType', 'approxllr','NoiseVariance',N0);
            demodulated_bits{u} = qamdemod(estimated_symbols{u}.',M,OutputType='bit');
            
            
            % LDPC decode the received bits
            [decoded_bits_ldpc{u},actualnumiter, finalparitychecks] = ldpcDecode(demodulatedLLR_ldpc, decodercfg, maxnumiter);
         
            % Calculate number of bit errors
            num_errors_ldpc = sum(original_bits{u} ~= decoded_bits_ldpc{u});
            num_errors = sum([original_bits{n}; zeros(2, 1)] ~= demodulated_bits{u});
            
            % Calculate BER
            BER_ldpc(u) = num_errors_ldpc / numBits;
            BER(u) = num_errors / (numBits+2);
            
            % Calculate BER for each user
            BER_users_ldpc(u, experiment) = BER_ldpc(u);
            BER_users(u, experiment) = BER(u);
        end     
    end
    average_BER_users_ldpc = mean(mean(BER_users_ldpc, 2));
    average_BER_users = mean(mean(BER_users, 2));

    BER_results_16QAM(1, l) = average_BER_users_ldpc;
    BER_results_16QAM(2, l) = average_BER_users;
end

