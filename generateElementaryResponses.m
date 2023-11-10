% PFC - Thiago Carvalho Bittencourt - EE - 2023 %

function aj = generateElementaryResponses(nFFT, fs, J)
    %initialize aj with zeros
    aj = zeros(J, nFFT);
    
    % Creates the frequency vector
    freq_vector = linspace(30, fs/2, nFFT/2 + 1);
    
    % Converting to Mel Scale
    mel_vector = 2595 * log10(1 + freq_vector / 700);
    
    % Uniformly spaced points in the Mel Frequency Scale
    mel_points = linspace(min(mel_vector), max(mel_vector), J + 2);
    
    % Convert the points back to Hz
    freq_points = 700 * (10.^(mel_points / 2595) - 1);
    
    % Convert the Frequency vector to the corresponding frequency bins
    bin_points = round(freq_points / fs * nFFT);
    
    % Loop for each triangular response
    for j = 1:J
       
        lower = bin_points(j);
        center = bin_points(j + 1);
        upper = bin_points(j + 2);
        
        lower = max(1, lower);
       
        rise = linspace(0, 1, center - lower + 1);
        
        fall = linspace(1, 0, upper - center + 1);
        
        aj(j, lower:upper) = [rise, fall(2:end)];
    end
end
