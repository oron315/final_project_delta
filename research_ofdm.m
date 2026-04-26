close all
clear
clc


% 1. Open the file
filename = 'raw_samples_12_Feb_2026_09_30_39_442_fs_10MHz.32fc';
fid = fopen(filename, 'r');

% 2. Read the data as 32-bit floats
% fread reads data into a 1D column vector
data = fread(fid, inf, 'float32');

% 3. Close the file
fclose(fid);

% 4. Reshape into Complex (Interleaved I/Q)
% Every pair is I + jQ
complex_data = data(1:2:end) + 1j*data(2:2:end);

% data inputs
fs = 10e6;
scs = 15e3;
T_sym = 1/scs;

M = 4;
nfft = 1024; % Number of FFT points
Ncp = 72; % Cyclic prefix length
Ncp_prime = 80;

% ----------- basic plots ------------ %
f = -fs/2:fs/length(complex_data):fs/2 - fs/length(complex_data);
t = 0:1/fs:length(complex_data)/fs-1/fs;

figure;
plot(t,complex_data)

figure;
plot(f,abs(fftshift(fft(complex_data))))

figure;
spectrogram(complex_data)

figure;
plot(real(complex_data),imag(complex_data))

% ---------------resample and stating to proccess the data---------------------%

data_upsample = resample(complex_data,1536,1000);
data_upsample = data_upsample.';
fs = nfft * scs;



% % ------------------ trying to padd with zeros --------------------- %
% data_upsample = [zeros(), fft(complex_data),zeros()]


% ---------create zadof chu------- %
Nzc = 601;

root = 600;
seq_4 = zadof_ofdm(Nzc,root);

root = 147;
seq_6 = zadof_ofdm(Nzc,root);



% -----------sync time--------------%

sync_time_4 = conv(data_upsample, flip(conj(seq_4)), "valid");

sync_time_6 = conv(data_upsample, flip(conj(seq_6)), "valid");

[~,max_sync_time_4] = max(sync_time_4);
[~,max_sync_time_6] = max(sync_time_6);

% sync_data = data_upsample(max_sync_time-3*(nfft+Ncp)-8:max_sync_time+6*(nfft+Ncp)+8);
sync_data = data_upsample(max_sync_time_4-3*(nfft+Ncp)-Ncp_prime:max_sync_time_4+6*nfft+4*Ncp+Ncp_prime);

t = 0:1/fs:length(sync_data)/fs-1/fs;

figure;
plot(abs(sync_time_4))
title('The correlation of zadoff chu 4')
figure;
plot(abs(sync_time_6))
title('The correlation of zadoff chu 6')

figure;
plot(t,abs(sync_data))

% -----------sync freq--------------%
% f = -fs/2:fs/nfft:fs/2 - fs/nfft;
% t = 0:1/fs:length(sync_data)/fs-1/fs;
% 
% 
% [~,sync_freq_idx] = max(abs(fftshift(fft(seq_4 .* conj(sync_data((1024+Ncp)*3+8+1:(1024+Ncp)*4+8-Ncp))))));
% [~,sync_freq_idx] = max(abs(fftshift(fft(seq_6 .* conj(sync_data((1024+Ncp)*5+8+1:(1024+Ncp)*6+8-Ncp))))));
% 
% sync_freq_big = f(sync_freq_idx); % big freq offset
% 
% % finding out that the freq offset is less than sub carrier spacing freq
% figure;
% plot(abs(fftshift(fft(seq_4 .* conj(sync_data((1024+Ncp)*3+8+1:(1024+Ncp)*4+8-Ncp))))))
% 
% figure;
% plot(abs(fftshift(fft(seq_6 .* conj(sync_data((1024+Ncp)*5+8+1:(1024+Ncp)*6+8-Ncp))))))

% -----------correcting soft freq-------------- %
% sync_freq_soft = angle(mean(seq_6 .* conj(sync_data((nfft+Ncp)*5+8+1:(nfft+Ncp)*6+8-Ncp)).*conj(seq_4 .* conj(sync_data((nfft+Ncp)*3+8+1:(1024+Ncp)*4+8-Ncp)))))*scs/(2*pi*2)


sync_freq_soft = angle(mean(seq_6 .* conj(sync_data((nfft+Ncp)*5+Ncp_prime+1:(nfft+Ncp)*5+Ncp_prime+nfft).*conj(seq_4 .* conj(sync_data((nfft+Ncp)*3+Ncp_prime+1:(nfft+Ncp)*3+nfft+Ncp_prime))))))*scs/(2*pi*2)
% sync_freq_soft_better = diff(angle(seq_6 .* conj(sync_data((nfft+Ncp)*5+8+1:(nfft+Ncp)*6+8-Ncp))))*fs/(2*pi);


sync_freq_data = sync_data .* exp(-1j*2*pi*(sync_freq_soft).*t);

% ---------demod ofdm to qpsk symbols---------- %

% ------ remove cp ----------%
sync_freq_data_first_ofdm = sync_freq_data(Ncp_prime+1:nfft+Ncp_prime);
sync_freq_data_last_ofdm = sync_freq_data(nfft*8+Ncp*9+8*2+1:nfft*9+Ncp*9+8*2);

sync_data_middle = reshape(sync_freq_data(nfft+Ncp_prime+1:nfft*8+Ncp*8+8),Ncp+nfft,[]);


sync_data_middle_mat = sync_data_middle(Ncp+1:nfft+Ncp,:);

sync_data_mat = [sync_freq_data_first_ofdm.',sync_data_middle_mat,sync_freq_data_last_ofdm.'];

% --------------------------- %

symbols_mat = fftshift(fft(sync_data_mat));

symbols = reshape(symbols_mat,1,[]);


figure;
scatter(real(symbols_mat(213:813,[2,3,5,7,9])),imag(symbols_mat(213:813,[2,3,5,7,9])))


figure;
plot(abs(symbols_mat(:,5)))



% ----- fix phase ------%






function seq_freq = zadof_ofdm(Nzc,root)
    
    m = 0:Nzc-1;
    seq = exp((-1j*pi*root*m.*(m+1)/Nzc));

   seq = [zeros(1,212),seq,zeros(1,1024-813)];
  
   seq_freq = ifft(ifftshift(seq));

end

