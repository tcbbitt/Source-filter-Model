% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
 function e_n_t = generateExcitationSpectrum(pitchData, numInstruments, fs, nFFT)

% Initialize the cell array "e_n_t" with zeros
 e_n_t = cell(1, size(pitchData, 2)); 

 for y=1:size(pitchData, 2)

     e_n_t{1,y} = zeros(numInstruments,nFFT); 

 end

    %Generate frequency vector
    f = linspace(0, fs / 2, nFFT); % from 0 to fs/2 with nFFT points
    nFrames = size(pitchData, 2);


    for t = 1:nFrames
        for n = 1:numInstruments
           % calculate the frequencies 
           harmonicFreqs = pitchData(n,t):pitchData(n,t):fs/2;  
           
           for h = 1:length(harmonicFreqs)
               % check the index of the closest frequency bin
              [~,cidx] = min(abs(f(2:end) - harmonicFreqs(h)));
               % set to 1
              e_n_t{1,t}(n,cidx) = 1;
           end
        end
    end
 end


