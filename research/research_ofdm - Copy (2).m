close all
clear
clc

%% ----------Open the file----------%
filename = 'raw_samples_12_Feb_2026_09_30_39_442_fs_10MHz.32fc';
fid = fopen(filename, 'r');

% Read the data as 32-bit floats
% fread reads data into a 1D column vector
data = fread(fid, inf, 'float32');

% Close the file
fclose(fid);

% Reshape into Complex (Interleaved I/Q)
% Every pair is I + jQ
complex_data = data(1:2:end) + 1j*data(2:2:end);

%% -----data inputs---------%
fs_start = 10e6;
scs = 15e3;
T_sym = 1/scs;

nfft = 1024; % Number of FFT points
Ncp = 72; % Cyclic prefix length
Ncp_prime = 80;

%% ----------- basic plots ------------ %
f = -fs_start/2:fs_start/length(complex_data):fs_start/2 - fs_start/length(complex_data);
t = 0:1/fs_start:length(complex_data)/fs_start-1/fs_start;

figure;
plot(t,complex_data)

figure;
plot(f,abs(fftshift(fft(complex_data))))

%% ---------------resample and stating to proccess the data---------------------%
% new sampling rate
fs = nfft * scs;

data_upsample = resample(complex_data,fs,fs_start);                  
data_upsample = data_upsample.';

%% ---------create zadof chu------- %
Nzc = 601;

root = 600;
seq_4 = zadof_ofdm(Nzc,root);

root = 147;
seq_6 = zadof_ofdm(Nzc,root);

%% -----------sync time--------------%

sync_time_4 = conv(data_upsample, flip(conj(seq_4)), "valid");

sync_time_6 = conv(data_upsample, flip(conj(seq_6)), "valid");


[~,max_sync_time_4] = max(sync_time_4);
[~,max_sync_time_6] = max(sync_time_6);

sync_data = data_upsample(max_sync_time_4-3*(nfft+Ncp)-Ncp_prime:max_sync_time_4+6*nfft+4*Ncp+Ncp_prime);

figure;
plot(abs(sync_time_4))
title('The correlation of zadoff chu 4')
figure;
plot(abs(sync_time_6))
title('The correlation of zadoff chu 6')

t = 0:1/fs:length(sync_data)/fs-1/fs;
figure;
plot(t,abs(sync_data))
title("The packet")

%% -----------sync freq low resolution--------------%
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

%% -----------correcting soft freq-------------- %

margin = 20;
sym_6 = sync_data((nfft+Ncp)*5+Ncp_prime+1-margin:(nfft+Ncp)*5+Ncp_prime+nfft+margin);

sync_freq_soft = angle(mean(seq_6 .* conj(sync_data((nfft+Ncp)*5+Ncp_prime+1:(nfft+Ncp)*5+Ncp_prime+nfft).*conj(seq_4 .* conj(sync_data((nfft+Ncp)*3+Ncp_prime+1:(nfft+Ncp)*3+nfft+Ncp_prime))))))*scs/(2*pi*2);
% sync_freq_soft_better = diff(angle(seq_6 .* conj(sync_data((nfft+Ncp)*5+8+1:(nfft+Ncp)*6+8-Ncp))))*fs/(2*pi);
% sync_freq_kay =  kay(seq_6,(sync_data((nfft+Ncp)*5+Ncp_prime+1:(nfft+Ncp)*5+Ncp_prime+nfft)),fs);

% sync_freq_data = sync_freq_data .* exp(-1j*2*pi*sync_freq_kay.*t);

% fixing freq error twice
freq_offset_cp = cp(sync_data,Ncp,Ncp_prime,fs)

sync_freq_data = sync_data .* exp(-1j*2*pi*(-freq_offset_cp).*t);


freq_offset_cp_2 = cp(sync_freq_data,Ncp,Ncp_prime,fs);

sync_freq_data = sync_freq_data .* exp(-1j*2*pi*(-freq_offset_cp_2).*t);

%% ---------demod ofdm to QPSK symbols---------- %

fixed_symbols_mat = ofdm_demod(sync_freq_data, Ncp,Ncp_prime,seq_6,fs);


%% ----- fix time once more ------%
t = 0:1/fs:length(data_upsample)/fs-1/fs;

data_fixed_freq = data_upsample .* exp(-1j*2*pi*(-freq_offset_cp).*t).* exp(-1j*2*pi*(-freq_offset_cp_2).*t);

sync_time_4 = conv(data_fixed_freq, flip(conj(seq_4)), "valid");

sync_time_6 = conv(data_fixed_freq, flip(conj(seq_6)), "valid");

[~,max_sync_time_40] = max(sync_time_4);
[~,max_sync_time_60] = max(sync_time_6);

sync_data = data_fixed_freq(max_sync_time_40-3*(nfft+Ncp)-Ncp_prime:max_sync_time_40+6*nfft+4*Ncp+Ncp_prime);

fixed_symbols_mat = ofdm_demod(sync_data, Ncp,Ncp_prime,seq_6,fs);

M = 4;
hex_msg = qpsk_demod(fixed_symbols_mat, M);



%% -----functions------%

function seq_freq = zadof_ofdm(Nzc,root)
    
    m = 0:(Nzc-1);
    seq = exp((-1j*pi*root*m.*(m+1)/Nzc));

    seq = [zeros(1,212),seq,zeros(1,(1024-813))];
  
    seq_freq = ifft(ifftshift(seq));

end


function freq_offset = kay(zadof,pilot,fs)
    L = length(pilot);
    k = (1:L-1);
    z = conj(zadof) .* pilot;

    z2 = z(1:end - 1);
    z1 = z(2:end);
    
    w = (3/2) * (L / (L^2 - 1)) * (1 - ((2 * k - L) / L).^2);

    freq_offset = sum(w .* (angle(z1 .* conj(z2)))) * fs / (2 * pi);
end


function freq_offset = cp(packet, Ncp,Ncp_prime,fs)
   
    nfft = 1024;
    sync_data_middle = reshape(packet(nfft+Ncp_prime+1:nfft*8+Ncp*7+Ncp_prime),Ncp+nfft,[]);

    cps_start = sync_data_middle(1:Ncp,:);

    cps_end = sync_data_middle(end-Ncp+1:end,:);

    match = cps_start .* conj(cps_end);

    figure;
    plot(angle(match))

    figure;
    plot(mean(diff(angle(match)),2))

    avg = angle(mean(match(20:50,:),"all"));

    freq_offset = avg*fs/(2*pi*1096);
end


function phase_offset = find_phase_offset(zadof,pilot,fs)
        
    % zadof meaning the real zadof chu and pilot is the symbol in the
    % packet with that zadof chu
    phase_offset = angle(mean(conj(zadof) .* pilot));
    
end

% 
% function H_est = estimate_channel(zadof,pilot)
% 
%     H = pilot .* zadof.';
%     H_est = movmean(H,64,'Endpoints','shrink');
% end



function find_freq_with_matrix(pilot_with_margin,zadof,fs)

    f_scan = (-1000:0.1:1000).';
    t = 0:1/fs:length(pilot_with_margin)/fs-1/fs;
    
    fft_mat = repmat(pilot_with_margin,length(f_scan),1) .* exp(-1j*2*pi*f_scan * t);

    cor_mat = conv2(fft_mat, zadof,"same");
    
    figure;
    % surf(t,f_scan.',abs(cor_mat));
    % hold on
    imagesc(t,f_scan.',abs(cor_mat).^12)
    colormap("jet");
   
    pow_mat = abs(cor_mat).^12;

    [~,idx] = max(pow_mat(:));
    [row,col] = ind2sub(size(pow_mat),idx);

    figure;
    contour(t,f_scan.',pow_mat)
    colormap
    colorbar
    title("contour map")

end


function fixed_symbols_mat = ofdm_demod(sync_freq_data, Ncp,Ncp_prime,zadof,fs)
    
    % take the sixth symbol from the packet
    nfft =1024;
    sym_6 = sync_freq_data((nfft+Ncp)*5+Ncp_prime+1:(nfft+Ncp)*5+Ncp_prime+nfft);

    % ------ forming the grid ----------%
    
    sync_freq_data_first_ofdm = sync_freq_data(Ncp_prime+1:nfft+Ncp_prime);
    sync_freq_data_last_ofdm = sync_freq_data(nfft*8+Ncp*9+8*2+1:nfft*9+Ncp*9+8*2);
    
    sync_data_middle = reshape(sync_freq_data(nfft+Ncp_prime+1:nfft*8+Ncp*8+8),Ncp+nfft,[]);
    
    % remove the cp
    sync_data_middle_mat = sync_data_middle(Ncp+1:nfft+Ncp,:);
    
    sync_data_mat = [sync_freq_data_first_ofdm.',sync_data_middle_mat,sync_freq_data_last_ofdm.'];
    
    symbols_mat = fftshift(fft(sync_data_mat),1);
  
    figure;
    scatter(real(symbols_mat(213:813,[2,3,5,6,7,9])),imag(symbols_mat(213:813,[2,3,5,6,7,9])))
    title("constellation QPSK")

    phase_offset = find_phase_offset(zadof,sym_6,fs);


    symbols_mat = symbols_mat * exp(-1j*phase_offset);

    figure;
    scatter(real(symbols_mat(213:813,[2,3,5,7,9])),imag(symbols_mat(213:813,[2,3,5,7,9])))

    figure;
    plot(abs(fftshift(fft((symbols_mat(:,3).^4).'))))

    
    %--------- fixing channel ---------- %

    % H_est = estimate_channel(symbols_mat(:,6),fftshift(fft(zadof)));


    h = (symbols_mat(:,6) ./ fftshift(fft(zadof)).');

    fixed_symbols_mat = symbols_mat ./ repmat(h,1,9);

    figure;
    scatter(real(fixed_symbols_mat(213:813,[2,3,5,7,8,9])),imag(fixed_symbols_mat(213:813,[2,3,5,7,8,9])))


end


function hex_msg = qpsk_demod(fixed_symbols_mat, M)

        % -----------QPSK demod-----------%

    fixed_symbols = fixed_symbols_mat([213:511,513:813],[2,3,5,7,8,9]);

    fixed_symbols = fixed_symbols(:);

    % Perform QPSK demodulation on the fixed symbols
    demodulated_symbols = pskdemod(conj(fixed_symbols),M,-pi/4,'gray',OutputType='bit');

    hex_msg = binaryVectorToHex(demodulated_symbols.');

end




