% Adding Gaussian noise to ECG signal
function noisy_ecg = add_gaussian_noise(Original_signal, snr_db)
signal_power = rms(Original_signal)^2;
noise_power = signal_power / (10^(snr_db/10));
noise = sqrt(noise_power) * randn(size(Original_signal));
noisy_ecg = Original_signal + noise;
end
 
% Load ECG signal data
load('/MATLAB Drive/100m.mat');
 
val = (val - 0) / 200; 
sig = val(1, 1:3600); 
 
Fs = 360; 
 
t = (0:length(sig)-1) / Fs; 
 
snr_db = 10; 
noisy_sig = add_gaussian_noise(sig, snr_db);
 
% Plot original and noisy ECG signals
figure;
subplot(2,1,1);
plot(t, sig, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ECG Amplitude (mV)');
title('Original ECG Signal');
grid on;
 
subplot(2,1,2);
plot(t, noisy_sig, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ECG Amplitude (mV)');
title('Noisy ECG Signal');
grid on;

 
% Variational Mode Decomposition (VMD)
alpha = 2000;      
tau = 0;            
K = 4;              
DC = 0;             
init = 1;           
tol = 1e-7;      

%Perform VMD
 
[u, ~] = VMD(noisy_sig, alpha, tau, K, DC, init, tol);

 

u = u';  
 

waveletName = 'db4'; 
for i = 1:K
    [c, l] = wavedec(u(:, i), 5, waveletName); 
    sigma = median(abs(c)) / 0.6745;
    thr = sigma * sqrt(2*log(length(c)));
    c = wthresh(c, 's', thr);

    % Reconstruction
    u(:, i) = waverec(c, l, waveletName);
end
 
% Reconstruct the signal from the denoised modes
denoised_sig = sum(u, 2);  
 
% Plot denoised ECG signal
figure;
plot(t, denoised_sig, 'LineWidth', 2);
xlabel('Time (s)');
ylabel('ECG Amplitude (mV)');
title('Denoised ECG Signal (VMD + Wavelet Denoising)');
grid on;

