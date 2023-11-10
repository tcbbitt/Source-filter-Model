% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
function cij = updateCij(r_t, g_n_i_t, e_n_t, cij, a_j)
   % Get Dimensions
    [I, N, T, K] = size(g_n_i_t);
    J = size(a_j, 1);
    L = size(cij{1, 1}, 2);
    
    % Update cij
    for i = 1:I
        for j = 1:J
            numerator = zeros(1, L);
            denominator = zeros(1, L);
            
            for n = 1:N
%                              
                for t = 1:T
%                     
                    for k = 1:K
%                         
                        % Update numerator and denominator 
                     
                        numerator = numerator + r_t(k, t) .* g_n_i_t(i, n, t, k) .* e_n_t{1,t}(n, k) .* a_j(j, k);
                        denominator = denominator + g_n_i_t(i, n, t, k) .* e_n_t{1,t}(n, k).* a_j(j, k);
                    end
                end
            end
            
            % Avoid division by zero and numerical error
            zeroIndices = (denominator == 0);
            denominator(zeroIndices) = 1e-10;
            
            % Update the cell cij{i, j}
            cij{i, j} = cij{i, j} * (numerator / denominator);
        end
    end
% guarantees that cij is non-negative
 for i = 1:I
        for j = 1:J
            cij{i,j} = abs(cij{i,j});
        end
 end
    
end
