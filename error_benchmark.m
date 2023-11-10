                %% error benchmark %%
% PFC - Thiago Carvalho Bittencourt - EE - 2023 %
clc;
clear;

%Select music
music = '09-Jesus';
%% Load data

load (['C:\Th\EE\0_PFC\reconstructed_files\' music '\var_' music '_int10_j24.mat']);

% define audio file path
audioFilePath = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '.wav']);
audioFilePath_1 = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '-violin.wav']);
audioFilePath_2 = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '-saxphone.wav']);
audioFilePath_3 = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '-clarinet.wav']);
audioFilePath_4 = (['C:\Th\EE\0_PFC\Bibliografia\2 - Database\Bach10\Bach10_v1.1\' music '\' music '-bassoon.wav']);

% number of test frames
testframes = 500;

%% Load and validate input mixture audio

% Original Music
% test file path
try
    [x_original, ~] = audioread(audioFilePath);
catch
    error('Could not read audio data from %s', audioFilePath);
end
x_original = x_original(:);

% Calculate the number of samples needed
numSamplesNeeded = 1 + (testframes - 1) * hopSizeSamples + windowLengthSamples;

% Keep only the first 'numSamplesNeeded' samples of the audio to test
x_original = x_original(1:numSamplesNeeded);

%% instrument 1
try
    [x_original_1, ~] = audioread(audioFilePath_1);
catch
    error('Could not read audio data from %s', audioFilePath_1);
end
x_original_1 = x_original_1(:); 

% Keep only the first 'numSamplesNeeded' samples of the audio to test
x_original_1 = x_original_1(1:numSamplesNeeded);

%% Instrument 2
try
    [x_original_2, ~] = audioread(audioFilePath_2);
catch
    error('Could not read audio data from %s', audioFilePath_2);
end
x_original_2 = x_original_2(:);  

% Keep only the first 'numSamplesNeeded' samples of the audio to test
x_original_2 = x_original_2(1:numSamplesNeeded);

%% Instrument 3
try
    [x_original_3, ~] = audioread(audioFilePath_3);
catch
    error('Could not read audio data from %s', audioFilePath_3);
end
x_original_3 = x_original_3(:);  % Make it a column vector

% Keep only the first 'numSamplesNeeded' samples of the audio to test
x_original_3 = x_original_3(1:numSamplesNeeded);

%% Instrument 4
try
    [x_original_4, ~] = audioread(audioFilePath_4);
catch
    error('Could not read audio data from %s', audioFilePath_4);
end
x_original_4 = x_original_4(:); 

% Keep only the first 'numSamplesNeeded' samples of the audio to test
x_original_4 = x_original_4(1:numSamplesNeeded);

% 
%% Original Music

%Calculate the noise signal
noise = x_original - reconstructed_normalized(1:numSamplesNeeded);

%Calculate the power of the noise
noisePower = sum(noise.^2);

%Calculate the power of the original signal
signalPower = sum(x_original.^2);

%Calculate uSDR
uSDR = 10 * log10(signalPower / noisePower);

%Print the average uSDR
fprintf('The average uSDR of the main music is: %f dB\n', uSDR);

%% instrument 1
%Calculate the noise signal
noise_1 = x_original_1 - reconstructed_normalized_instruments{1}(1:numSamplesNeeded);

%Calculate the power of the noise
noisePower_1 = sum(noise_1.^2);

%Calculate the power of the original signal
signalPower_1 = sum(x_original_1.^2);

%Calculate uSDR
uSDR_1 = 10 * log10(signalPower_1 / noisePower_1);

%Print the average uSDR
fprintf('The average uSDR of instrument 1 is: %f dB\n', uSDR_1);

%% instrument 2
%Calculate the noise signal
noise_2 = x_original_2 - reconstructed_normalized_instruments{2}(1:numSamplesNeeded);

%Calculate the power of the noise
noisePower_2 = sum(noise_2.^2);

%Calculate the power of the original signal
signalPower_2 = sum(x_original_2.^2);

%Calculate uSDR
uSDR_2 = 10 * log10(signalPower_2 / noisePower_2);

%Print the average uSDR
fprintf('The average uSDR of instrument 2 is: %f dB\n', uSDR_2);
% 
%% instrument 3
%Calculate the noise signal
noise_3 = x_original_3 - reconstructed_normalized_instruments{3}(1:numSamplesNeeded);

%Calculate the power of the noise
noisePower_3 = sum(noise_3.^2);

% Calculate the power of the original signal
signalPower_3 = sum(x_original_3.^2);

% Calculate SNR
uSDR_3 = 10 * log10(signalPower_3 / noisePower_3);

% Print the average SNR
fprintf('The average uSDR of instrument 3 is: %f dB\n', uSDR_3);

%% instrument 4
%Calculate the noise signal
noise_4 = x_original_4 - reconstructed_normalized_instruments{4}(1:numSamplesNeeded);

%Calculate the power of the noise
noisePower_4 = sum(noise_4.^2);

%Calculate the power of the original signal
signalPower_4 = sum(x_original_4.^2);
% 
%Calculate SNR
uSDR_4 = 10 * log10(signalPower_4 / noisePower_4);

%Print the average SNR
fprintf('The average uSDR of instrument 4 is: %f dB\n', uSDR_4);

%% Plot 
% create time vector
time = (0:(length(x_original)-1)) / fs;
% plot selected instrument over reference signal
plot(time,reconstructed_normalized_instruments{2}(1:numSamplesNeeded),'b');
xlabel('Tempo (s)');
ylabel('Amplitude Normalizada');
hold on;
plot(time, x_original_2, 'r')
hold on;
legend('Sinal reconstruído','Sinal original');
title('Sinal')
