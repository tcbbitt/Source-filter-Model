% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
function Gnit = updateG_n_i_t(r_t, e_n_t, cij, Gnit, aj)
    % Get Dimensions
    [I, N, T, K] = size(Gnit);
    [I, J] = size(cij);
    % Update Gnit
    for i = 1:I
        for n = 1:N
            for t = 1:T
                for k =1:min(K, size(e_n_t{t}, 2))
                    numerator = 0;
                    denominator = 0;
                    
                    if k > size(e_n_t{t}, 2)
                        warning('Index k exceeds array bounds for e_n_t{t}.');
                        continue;  
                    end
                   
                    for j = 1:J
                         %debug check
                        % Check bounds before the calculation
                        if k > size(aj, 2) || n > size(e_n_t{t}, 2) || t > size(r_t, 1)
                            warning(['out of bounds']);
                            continue;
                        end
                        numerator = numerator + r_t(k, t) .* e_n_t{1,t}(n, k) .* cij{i,j}(1, k) .* aj(j, k);
                        denominator = denominator + e_n_t{1,t}(n, k) .* cij{i,j}(1, k) .* aj(j, k);
                    end
            % avoid numerical error
            zeroIndices = (denominator == 0);
            denominator(zeroIndices) = 1e-10;

                    Gnit(i, n, t, k) = Gnit(i, n, t, k) * (numerator / denominator);
                  
                end
            end
        end
    end
    %guarantees non-negative values
    Gnit = abs(Gnit);
  end
