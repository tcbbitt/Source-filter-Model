% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
%min-max normalization
function normalized_data = min_max_normalization(data, min_val, max_val, new_min, new_max)
  
    % Perform Min-Max Normalization
    normalized_data = ((data - min_val) / (max_val - min_val)) * (new_max - new_min) + new_min;
end
