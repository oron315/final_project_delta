close all
clear
clc

%% testing SDR

% fs_1spc = 15.36e6;
% BW = fs_1spc * 2;

Ncp = [72 80];
nfft = 1024;

% time known data
packet_time = 643e-6;
time_between_packets = 640e-3;

fc_list = 1e9 * [2.3995, 2.4145,2.4295,2.4445,2.4595, 5.7965,5.7765, 5.7565];


%% ---------create zadof chu------- %
Nzc = 601;


root = 600;
seq_4 = zadof_ofdm(Nzc,root);

root = 147;
seq_6 = zadof_ofdm(Nzc,root);

%% ---------- system test 1 channel only----------- %

fc = fc_list(3);
fs = 15.36e6;

% fc = 98e6;
% fs = 15.36e6;
% taking 2 packets in time
time2packets = 5*(packet_time * 2 + time_between_packets);

% -------- sampling using SDR --------- %

% packet length in samples for fs=15.36e6
N_packet = (nfft * 9 + Ncp(1) * 7 + Ncp(2) * 2);

% resetting to packet not found
is_packet = false;

% keep sampling until an entire packet is detcted
while ~is_packet
    [is_packet,data,max_idx_sync_time_4] = process_1_channel(fc,fs,time2packets,N_packet,nfft,Ncp,seq_4);
end

% [is_packet,data,max_idx_sync_time_4,rx,scope] = process_1_channel(fc,fs,time2packets,N_packet,nfft,Ncp,seq_4);
% taking only the packet from the data
sync_data = data(max_idx_sync_time_4-3*(nfft+Ncp(1))-Ncp(2):max_idx_sync_time_4+6*nfft+4*Ncp(1)+Ncp(2));

% fix freq offset twice cause good
[fixed_packet,~] = cp_fix_freq(sync_data, Ncp(1),Ncp(2),fs);

[fixed_packet,~] = cp_fix_freq(fixed_packet, Ncp(1),Ncp(2),fs);

% Demod the ofdm symbols
fixed_symbols_mat = ofdm_demod(sync_freq_data, Ncp,Ncp_prime,seq_6,fs);

% Demod QPSK
M = 4;
bin_msg = qpsk_demod(fixed_symbols_mat, M);

%% process only one channel

function [is_packet,data,max_idx_sync_time_4] = process_1_channel(fc,fs,time2packets,N_packet,nfft,Ncp,seq_4)
    % gets data of length 2*packet_length 
    % returns boolean of there is a packet there and it is not to late or
    % early

    [data] = capture_samples(fc,fs,time2packets);
    
    [found_pilot,max_idx_sync_time_4] = find_pilot(data,seq_4);

    % borders for getting an entire packet of info
    min_idx_4 = 3 * nfft + 3 * Ncp(1) + Ncp(2);
    max_idx_4 = 2 * N_packet - (nfft * 6 + Ncp(1) * 5 + Ncp(2));
    
    is_packet = false;
    if found_pilot == 1 && ~((max_idx_sync_time_4 > max_idx_4) || (max_idx_sync_time_4 < min_idx_4))
        is_packet = true;
    end    
end

%% Finding pilot

function [found_pilot,max_idx_sync_time_4] = find_pilot(data,seq_4)
    
    found_pilot = 0;
    sync_time_4 = conv(data, flip(conj(seq_4)), "valid");
    
    [max_corr_4,max_idx_sync_time_4] = max(sync_time_4);

    % Calculate the peak to average ratio (PAR)
    par = (max_corr_4 / mean(abs(data)));

    if abs(par) > length(seq_4)
        found_pilot = 1;
    end
    
end

%% Sampling using the SDR

function [data1] = capture_samples(fc,fs,measure_time)

    rx = comm.SDRuReceiver(...
                  Platform ="B210", ...
                  SerialNum ="3591273", ...
                  CenterFrequency =fc, ...
                  MasterClockRate =fs, ...
                  DecimationFactor =1);
    
    
    rx.ReceiveAntennaPort = 'TX/RX';
    rx.Gain = 20;
    [data,metadata]= capture(rx,measure_time,"Seconds");
    
    data = double(data);
    
    data1 = data ./(2^15);
    
    spectrogram(data1);
    % 
    % sampleRate = rx.MasterClockRate/rx.DecimationFactor;
    % scope = spectrumAnalyzer(SampleRate=sampleRate);
    % scope(data1)
    % release(scope);

end

%% creating zadof chu as an OFDM symbol
function seq_freq = zadof_ofdm(Nzc,root)
    
    m = 0:(Nzc-1);
    seq = exp((-1j*pi*root*m.*(m+1)/Nzc));

    seq = [zeros(1,212),seq,zeros(1,(1024-813))];
  
    seq_freq = ifft(ifftshift(seq));

    % seq_freq = resample(seq_freq,scale,1);

end

%% correcting frequency offset

function [fixed_packet,freq_offset] = cp_fix_freq(packet, Ncp,Ncp_prime,fs)
   
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

    t = 0:1/fs:length(packet)/fs-1/fs; 

    fixed_packet = packet .* exp(1j*2*pi*freq_offset.*t);
end

%% Demod OFDM

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

    h = (symbols_mat(:,6) ./ fftshift(fft(zadof)).');

    fixed_symbols_mat = symbols_mat ./ repmat(h,1,9);

    figure;
    scatter(real(fixed_symbols_mat(213:813,[2,3,5,7,8,9])),imag(fixed_symbols_mat(213:813,[2,3,5,7,8,9])))


end

%% Demod QPSK

function demodulated_symbols = qpsk_demod(fixed_symbols_mat, M)

    fixed_symbols = fixed_symbols_mat([213:512,514:813],[2,3,5,7,8,9]);

    fixed_symbols = fixed_symbols(:);

    % Perform QPSK demodulation on the fixed symbols
    demodulated_symbols = pskdemod(conj(fixed_symbols),M,-pi/4,'gray',OutputType='bit');

    % hex_msg = binaryVectorToHex(demodulated_symbols.')

end

%% use later - trying to sample many channel at the same time

% function [fc] = recieve_wideband_signal(data,BW,fs_1spc,fc_now, fc_list)
% 
%      % Calculate the correlation peak to average data ratio (PAR)
%     par = find_is_signal(data,BW,fs_1spc)
% 
%     if par > length(seq_4)
% 
%         find_bits(data)
% 
%     else
% 
%         % moving to next fc option
%         fc_next = fc_list(find(fc_list == fc_now) + 1);
% 
%     end
% 
% end



%% viewing real time spectrum analyzer
% rx = comm.SDRuReceiver(...
%               Platform ="B210", ...
%               SerialNum ="3591273", ...
%               CenterFrequency =fc, ...
%               MasterClockRate =fs, ...
%               DecimationFactor =1);
% 
% 
% rx.ReceiveAntennaPort = 'TX/RX';
% rx.Gain = 50;
% 
% specAnalyzer = spectrumAnalyzer('SampleRate',fs);
% 
% for i=1:100000
%     [data,len] = rx();
%     % [data,metadata]= capture(rx,1,"Seconds");
%     specAnalyzer(data);
% end
% 
% release(rx);
% release(specAnalyzer);