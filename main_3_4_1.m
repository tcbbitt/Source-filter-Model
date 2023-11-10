%% Instrument reconstruction using "source filter" method %%
    % PFC - Thiago Carvalho Bittencourt - EE - 2023 %
clc;
clear;
%Select music
music = 'BachChorale-2-01-1';

% Define file paths for pitch and audio data
pitchFilePath = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '-GTF0s.mat']);
audioFilePath = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '.wav']);

% Load and validate pitch data
try
    pitchData = load(pitchFilePath);
    pitchData = cell2mat(struct2cell(pitchData));
catch
    error('Could not load %s', pitchFilePath);
end

% Convert pitch data from MIDI to Hz using MIDI Toolbox
[numInstruments, numFrames] = size(pitchData);
for frameIdx = 1:numFrames
    for instrIdx = 1:numInstruments
        pitchData(instrIdx, frameIdx) = midi2hz(pitchData(instrIdx, frameIdx));
    end
end

% number of test frames
testframes = 5;

% Parameters for the iterative algorithm (NMF)
maxIterations = 100;
Iteration_c = 0;
 
% Keep only the first X frames of the pitch data
pitchData = pitchData(:, 1:testframes);
[numInstruments, numFrames] = size(pitchData);

%% Parameters
% Sampling and window parameters 
fs = 44100;                 % Sample rate in Hz
windowLengthSec = 0.046;    % Window length in seconds
hopSizeSec = 0.01;          % Hop size in seconds
numFilters = 24;            % Number of filters

% Convert time values to samples
windowLengthSamples = ceil(windowLengthSec * fs) + 1;
hopSizeSamples = round(hopSizeSec * fs);

% Create Hamming window
hammWindow = hamming(windowLengthSamples);

%% Load and validate input mixture audio
try
    [x, ~] = audioread(audioFilePath);
catch
    error('Could not read audio data from %s', audioFilePath);
end
x = x(:);

% Calculate the number of samples needed for the first testframes windows 
numSamplesNeeded = 1 + (testframes - 1) * hopSizeSamples + windowLengthSamples;

% Keep only the first 'numSamplesNeeded' samples of the audio to test
x = x(1:numSamplesNeeded);
%% Perform STFT
numWindows = floor((length(x) - windowLengthSamples) / hopSizeSamples) + 1;
stftResult = zeros(windowLengthSamples, numWindows);

for windowIdx = 1:numWindows
    startIdx = (windowIdx - 1) * hopSizeSamples + 1;
    endIdx = startIdx + windowLengthSamples - 1;
    
    % Validate indices
    if endIdx > length(x)
        break;
    end
    
    % Apply Hamming window and compute FFT
    windowedSignal = x(startIdx:endIdx) .* hammWindow;
    stftResult(:, windowIdx) = fft(windowedSignal);
end
%% Initialize variables for the algorithm

% Get the size of the FFT
nFFT = size(stftResult, 1);

% Get phase of the STFT
originalPhase = angle(stftResult);

% Initialize Cij as a Cell array with absolute values of Gaussian noise
Cij = initializeCij(numInstruments, numFilters, nFFT);

% Initialize g_n_i_t with absolute values of Gaussian noise
Gnit = abs(randn(numInstruments, numInstruments,numWindows, nFFT));

% Initialize reconstructed ratio
rt = abs(randn(size(stftResult)));

% Generate elementary responses (aj)

aj = generateElementaryResponses(nFFT, fs, numFilters);

% Calculate e_n_t for each frame and each instrument
e_n_t = generateExcitationSpectrum(pitchData, numInstruments, fs, nFFT);

% Initialize x_hat_i_t for each instrument
x_hat_i_t = cell(numInstruments,1);

%% Start the iterative algorithm
for iter = 1:maxIterations

    % Update Cij
     Cij = updateCij(rt, Gnit, e_n_t,Cij, aj); 
    
    % Update G_n_i_t
     Gnit = updateG_n_i_t(rt, e_n_t, Cij,Gnit, aj);
    
    % Calculate reconstructed signal 

     x_hat_t = calculate_x_hat_t(Cij, Gnit, e_n_t, aj); 

    % update rt
     for t = 1:numFrames
        rt(:, t) = abs(stftResult(:, t)) ./ abs(x_hat_t(:, t));  
     end
     %print the number of iterations for control
     Iteration_c = Iteration_c + 1;
     fprintf('number of iterations: %f', Iteration_c);
end

%% reconstruct x_hat for each instrument
 
     for idx = 1:numInstruments
     x_hat_i_t{idx} = calculate_x_hat_i_t(Cij, Gnit, e_n_t, aj, idx);    
     end

%% Perform inverse STFT to get time-domain signals 

% Initialize variables to store time-domain signals
%whole signal 
x_reconstructed = zeros(size(x));

%each instrument using yi,t(k)
 x_reconstructed_instruments = cell(numInstruments,1);
for i = 1:numInstruments
    x_reconstructed_instruments{i} = zeros(size(x));
end
y_i_t = cell(numInstruments, 1);
for i = 1:numInstruments
    y_i_t{i} = zeros(size(stftResult));
end

%each instrument using x_hat,i,t(k)
x_reconstructed_instruments_hat = cell(numInstruments,1);
for i = 1:numInstruments
    x_reconstructed_instruments_hat{i} = zeros(size(x));
end

% Pre allocate for speed
%Whole signal

 ifft_result = zeros(size(stftResult));
% cell array of instruments
 
 ifft_result_i = cell(1,numInstruments);

for i = 1:numInstruments
   ifft_result_i{i} = zeros(size(stftResult));
end

 ifft_result_i_hat = cell(1,numInstruments);

for i = 1:numInstruments
   ifft_result_i_hat{i} = zeros(size(stftResult));
end

%% Reconstruction

 
 for t=1:numFrames
     
     % combining with the original phase for inverse STFT (for x_hat,t)
        ifft_result(:, t) = ifft(abs(x_hat_t(:, t)) .* exp(1i * originalPhase(:, t)), 'symmetric');
        
        %time-domain signal
        startIdx = (t - 1) * hopSizeSamples + 1;
        endIdx = startIdx + windowLengthSamples - 1;
        x_reconstructed(startIdx:endIdx) = x_reconstructed(startIdx:endIdx) + ifft_result(:, t) .* hammWindow;
        x_reconstructed(startIdx:endIdx) = x_reconstructed(startIdx:endIdx) + ifft_result(:, t) .* hammWindow;
 
        % For each instrument
    for i = 1:numInstruments   
       
        % Compute y_{i,t}(k)
      
          y_i_t{i}(:, t) = (x_hat_i_t{i}(:, t) ./ x_hat_t(:, t)) .* abs(stftResult(:, t));
        
        % combining with the original phase for inverse STFT (for y,i,t and x_hat,i,t)
        %using y,i
        ifft_result_i{i}(:, t) = ifft(abs(y_i_t{i}(:, t)) .* exp(1i * originalPhase(:, t)), 'symmetric');
        %using x_hat_i
        ifft_result_i_hat{i}(:, t) = ifft(abs(x_hat_i_t{i}(:, t)) .* exp(1i * originalPhase(:, t)), 'symmetric');


        %time-domain signal
        startIdx = (t - 1) * hopSizeSamples + 1;
        endIdx = startIdx + windowLengthSamples - 1;
        %using y,i
        x_reconstructed_instruments{i}(startIdx:endIdx) = x_reconstructed_instruments{i}(startIdx:endIdx) + ifft_result_i{i}(:, t) .* hammWindow;
        %using x_hat,i
        x_reconstructed_instruments_hat{i}(startIdx:endIdx) = x_reconstructed_instruments_hat{i}(startIdx:endIdx) + ifft_result_i_hat{i}(:, t) .* hammWindow;
    end
 end
    
    
%% Normalization 
% Save as .wav
    
%For the whole signal
reconstructed_normalized = x_reconstructed / max(abs(x_reconstructed(:)));    
audiowrite(['C:\Th\EE\0_PFC\reconstructed_files\' music '\whole_signal_reconstructed_normalized.wav'], reconstructed_normalized, fs);
    
% For each instrument
for i = 1:numInstruments
    % Normalize the audio data
    %y,i,t(k)
    reconstructed_normalized_instruments{i,1} = x_reconstructed_instruments{i} / max(abs(x_reconstructed_instruments{i}(:)));
    %x_hat,i,t(k)
    reconstructed_normalized_instruments_hat{i,1} = x_reconstructed_instruments_hat{i} / max(abs(x_reconstructed_instruments_hat{i}(:)));
    
    % Save the normalized audio data
    %%y,i,t(k)
    audiowrite(['C:\Th\EE\0_PFC\reconstructed_files\' music '\instrument_' num2str(i) '_reconstructed_normalized.wav'], reconstructed_normalized_instruments{i,1}, fs);
     %%x+,i,t(k)
    audiowrite(['C:\Th\EE\0_PFC\reconstructed_files\' music '\instrument_' num2str(i) '_reconstructed_normalized_xhat.wav'], reconstructed_normalized_instruments_hat{i}, fs);
end

% %% Non normalized signal (if needed)
% % Save as .wav
% 
% %For the whole signal  
% audiowrite('C:\Th\EE\0_PFC\reconstructed_files\BachChorale-2-02-5\whole_signal_reconstructed.wav', x_reconstructed, fs);
% 
% % For each instrument
% for i = 1:numInstruments
% 
%     % Save the audio data
%     %%y,i,t(k)
%     audiowrite(['C:\Th\EE\0_PFC\reconstructed_files\BachChorale-2-02-5\instrument_' num2str(i) '_reconstructed.wav'], x_reconstructed_instruments{i,1}, fs);
%      %%x,i,t(k)
%     audiowrite(['C:\Th\EE\0_PFC\reconstructed_files\BachChorale-2-02-5\instrument_' num2str(i) '_reconstructed_xhat.wav'], x_reconstructed_instruments_hat{i}, fs);
% end

%% Save variables for further testing
save(['C:\Th\EE\0_PFC\reconstructed_files\' music '\var_' music '_int' num2str(maxIterations) '_j' num2str(numFilters)]);


