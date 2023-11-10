% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
function x_hat_i_t = calculate_x_hat_i_t(Cij, Gnit, e_n_t, aj, instrument_index)
    % Get dimensions
    [I, N, T, K] = size(Gnit);  % Instruments, Number of notes, Frames, FFT Length
    J = size(aj, 1);            % Number of filters

    % Initialize x_hat_i_t with zeros
    x_hat_i_t = zeros(T, K);

    % Specify instrument
    i = instrument_index;

    % Iterate through T, N, and K to generate x_hat_i_t
    for n = 1:N
        for k = 1:K
            for t = 1:T
                % Sum over J
                sum_j = 0;
                for j = 1:J
                    sum_j = sum_j + Cij{i, j}(1, k).* aj(j, k);
                end

                % Update x_hat_i_t for each k
                x_hat_i_t(t, k) = x_hat_i_t(t, k) + Gnit(i, n, t, k) .* e_n_t{1, t}(n, k) .* sum_j;

                % Prevent numerical errors
                if x_hat_i_t(t, k) == 0
                    x_hat_i_t(t, k) = 1e-10;
                end
            end
        end
    end
    %transpose to the right format
     x_hat_i_t = x_hat_i_t';
end
