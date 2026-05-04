close all
clear
clc

%% ------INPUT DATA------%
% time known data
packet_time = 643e-6;
time_between_packets = 640e-3;

Ncp = [72 80];
nfft = 1024;
fc_list = 1e9 * [2.3995, 2.4145,2.4295,2.4445,2.4595];

%% ---------create zadof chu------- %
Nzc = 601;

root = 600;
seq_4 = zadof_ofdm(Nzc,root);

root = 147;
seq_6 = zadof_ofdm(Nzc,root);

%% ---------- system test 1 channel only----------- %
% only one channel
fc = fc_list(3);
fs = 15.36e6;

% taking 5 packets in time
time5packets = (packet_time * 5 + time_between_packets * 4);

% -------- sampling using SDR --------- %

% packet length in samples for fs=15.36e6
N_packet = (nfft * 9 + Ncp(1) * 7 + Ncp(2) * 2);

% resetting to packet not found
is_packet = false;

% keep sampling until an entire packet is detcted
fc_next = fc_list(3);
while ~is_packet
    [is_packet,data,max_idx_sync_time_4] = process_1_channel(fc_next,fs,time_between_packets,nfft,Ncp,seq_4);
    % frequency hopping because there are 5-20 packets in each freq
    [fc_next] = hopping_between_freq(is_packet,fc_next, fc_list);
    message = "i just hopped! :)"
end

% taking only the packet from the data
% each row of sync_data is a different packet

% sync_data = zeros(length(max_idx_sync_time_4),N_packet+1);
% for i = 1:length(max_idx_sync_time_4)
%     sync_data(i,:) = data(max_idx_sync_time_4(i)-3*(nfft+Ncp(1))-Ncp(2):max_idx_sync_time_4(i)+6*nfft+4*Ncp(1)+Ncp(2));
% end

% packet
sync_data = data(max_idx_sync_time_4-3*(nfft+Ncp(1))-Ncp(2):max_idx_sync_time_4+6*nfft+4*Ncp(1)+Ncp(2));

% fix freq offset twice cause good
% [fixed_packet,~] = cp_fix_freq(sync_data(end,:), Ncp(1),Ncp(2),fs);

[fixed_packet,~] = cp_fix_freq(sync_data, Ncp(1),Ncp(2),fs);

% Demod the ofdm symbols
fixed_symbols_mat = ofdm_demod(fixed_packet, Ncp(1),Ncp(2),seq_6,fs);

% Demod QPSK
M = 4;
raw_bin = qpsk_demod(fixed_symbols_mat, M);
raw_bin = raw_bin.'; % must be row vector !!!

%% ---------- Amittai bits -------------%
r_xor_mask = xor_mask(raw_bin);
r_cyclic_buffer = cyclic_buffer(r_xor_mask);
r_turbo_and_intereaver = turbo_and_interleaver(r_cyclic_buffer);
r_crc = crc(r_turbo_and_intereaver);

proccesed_bin = r_crc(1:91*8);
proccesed_text = binaryVectorToHex(proccesed_bin);

keys = ["payload_length", "unknown1", "version", "sequence_number",...
    "states_info", "serial", "long", "lat", "altitude", "height",...
    "v_north", "v_east", "v_up", "unknown2", "time", "app_lat"...
    "app_long","home_long","home_lat", "device_type","uuid_len", "uuid"...
    "crc", "crc_valid"];

% a repetition at the end because crc_valid is based on crc
lengths = [1 1 1 2 2 16 4 4 2 2 2 2 2 2 8 4 4 4 4 1 1 20 2 2];
offsets = [0 1 2 3 5 7 23 27 31 33 35 37 39 41 43 51 55 59 63 67 68 69 89 89];

d = dictionary();
for i = (1:length(offsets))
    current_bits = proccesed_bin(8*offsets(i) + 1:8*offsets(i) + 8*lengths(i));
    % little endian
    current_bits = reshape(flip(reshape(current_bits, 8, []), 2), 1, []);
    switch i
        % this is wrong because some are unsighned and some not
        case {1 2 3 14 21}  % int 8
            value = bin2dec(num2str(current_bits, '%d'));
            if value > 32767
                value = value - 65536;
            end
            value = string(value);
        case {4 9 10 11 12 13} % int 16
            value = bin2dec(num2str(current_bits, '%d'));
            if value > 32767
                value = value - 65536;
            end
            value = string(value);
        case {7 8 16 17 18 19}  % degrees, int 32
            value = bin2dec(num2str(current_bits, '%d')) / 174533;
            if value > 32767
                value = value - 65536;
            end
            value = string(value);
        case 15  % time, int64
            value = bin2dec(num2str(current_bits, '%d'));
            value = string(datetime(value/1000, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC'));
        case {6 22}  % ascii, big endian
            current_bits = reshape(flip(reshape(current_bits, 8, []), 2), 1, []);
            value = string(char(bin2dec(reshape(char(current_bits + '0'), 8, []).'))'); 
        case 5  % flags
            value = strjoin(string(logical(current_bits)));
        case 20
            value = string(bin2dec(num2str(current_bits, '%d')));
        case {23 24}  % crc
            [~, err] = nrCRCDecode(proccesed_bin, "16");
            % THIS IS WRONG!
            if i == 23
                value = binaryVectorToHex(current_bits);
            else
                % value = string(sum(err) ~= 0);
                value = "true";
            end
            
        otherwise
            value = join(string(-1));
    end
    d = insert(d, keys(i), value);
    
end

% print
d

%% GUI MAP
geoplot(str2double(d("app_lat")), str2double(d("app_long")), 'r*', 'MarkerSize', 10)

% end of code
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process only one channel

function [is_packet,data,max_idx_sync_time_4] = process_1_channel(fc,fs,time,nfft,Ncp,seq_4)
    % gets data of length 2*packet_length 
    % returns boolean of there is a packet there and it is not to late or
    % early

    [data] = capture_samples(fc,fs,time);
    data = data.';
    % data = load('data_2.4295fc_the_best.mat').data;
    % data = data.';
    % 
    [is_packet,max_idx_sync_time_4] = find_packet(data,seq_4,nfft,Ncp);

end

%% Finding pilot

function [is_packet,max_corr_idx] = find_packet(data,seq_4,nfft,Ncp)
    % returns boolean if there is packet in the data and the start index
    
    % max correlation
    corr_zc4 = conv(data, flip(conj(seq_4)), "valid");
    [max_corr_val,max_corr_idx] = maxk(corr_zc4,1);   
    
    % PAPR threshold
    par = abs((max_corr_val ./ mean(abs(data))));
    threshold = length(seq_4);
    threshold = threshold / 2;

    % threshold is (max_corr_val > mean(abs(data)) * length(seq_4)) / 2
    found_pilot = par > threshold;

    % borders for getting an entire packet of info
    min_idx_4 = (3 * nfft + 3 * Ncp(1) + Ncp(2));
    max_idx_4 = (length(data) - (nfft * 6 + Ncp(1) * 5 + Ncp(2)));

    enough_packet_before = max_corr_idx < max_idx_4;
    enough_packet_after = max_corr_idx > min_idx_4;
    enough_packet = and(enough_packet_before,enough_packet_after);

    found_pilot_full = and(found_pilot,enough_packet);
    is_packet = sum(found_pilot_full) > 0;

    max_corr_idx = max_corr_idx(found_pilot_full==1); % find the idx the work
        
end

%% Sampling using the SDR

function [data1] = capture_samples(fc,fs,measure_time)

    rx = comm.SDRuReceiver(...
                  Platform ="B210", ...
                  SerialNum ="3557010", ...
                  CenterFrequency =fc, ...
                  MasterClockRate =fs, ...
                  DecimationFactor =1);
    
    
    rx.ReceiveAntennaPort = 'TX/RX';
    rx.Gain = 20;

    [data,~]= capture(rx,measure_time,"Seconds");
    
    data = double(data);
    data1 = data ./(2^15); 
end

%% creating zadof chu as an OFDM symbol
function seq_freq = zadof_ofdm(Nzc,root)
    
    m = 0:(Nzc-1);
    seq = exp((-1j*pi*root*m.*(m+1)/Nzc));

    seq = [zeros(1,212),seq,zeros(1,(1024-813))];
  
    seq_freq = ifft(ifftshift(seq));

end

%% correcting frequency offset
function [fixed_packet,freq_offset] = cp_fix_freq(packet, Ncp,Ncp_prime,fs)
   
    nfft = 1024;
    sync_data_middle = reshape(packet(nfft+Ncp_prime+1:nfft*8+Ncp*7+Ncp_prime),Ncp+nfft,[]);

    cps_start = sync_data_middle(1:Ncp,:);

    cps_end = sync_data_middle(end-Ncp+1:end,:);

    match = cps_start .* conj(cps_end);

    % figure;
    % plot(angle(match))
    % 
    % figure;
    % plot(mean(diff(angle(match)),2))

    avg = angle(mean(match(20:50,:),"all"));

    freq_offset = avg*fs/(2*pi*1024);

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
  
    % figure;
    % scatter(real(symbols_mat(213:813,[2,3,5,6,7,9])),imag(symbols_mat(213:813,[2,3,5,6,7,9])))
    % title("constellation QPSK")

    phase_offset = find_phase_offset(zadof,sym_6,fs);


    symbols_mat = symbols_mat * exp(-1j*phase_offset);

    % figure;
    % scatter(real(symbols_mat(213:813,[2,3,5,7,9])),imag(symbols_mat(213:813,[2,3,5,7,9])))
    % 
    % 
    %--------- fixing channel ---------- %

    h = (symbols_mat(:,6) ./ fftshift(fft(zadof)).');

    fixed_symbols_mat = symbols_mat ./ repmat(h,1,9);

    % figure;
    % scatter(real(fixed_symbols_mat(213:813,[2,3,5,7,8,9])),imag(fixed_symbols_mat(213:813,[2,3,5,7,8,9])))


end


function phase_offset = find_phase_offset(zadof,pilot,fs)
        
    % zadof meaning the real zadof chu and pilot is the symbol in the
    % packet with that zadof chu
    phase_offset = angle(mean(conj(zadof) .* pilot));
    
end

%% Demod QPSK

function demodulated_symbols = qpsk_demod(fixed_symbols_mat, M)

    fixed_symbols = fixed_symbols_mat([213:512,514:813],[2,3,5,7,8,9]);

    fixed_symbols = fixed_symbols(:);

    % Perform QPSK demodulation on the fixed symbols
    demodulated_symbols = pskdemod(conj(fixed_symbols),M,-pi/4,'gray',OutputType='bit');

    % hex_msg = binaryVectorToHex(demodulated_symbols.')

end


%% BITS FUNCTIONS
function interleaved = apply_p2_interleaver(data)
    % data length is 1408
    % interleavered length is 1408
    p2_interleaver = mod((43 * (1:1408) + 88 * (1:1408).^2), 1408);
    interleaved = zeros(1, length(data));
    for i = (1:length(p2_interleaver))
        interleaved(i) = data(p2_interleaver(i) + 1);
    end
end

function r_p2_interleaver = p2_interleaver(data)
    % data length is 1408
    % r_p2_interleaver length is 1408
    p2_interleaver = mod((43 * (1:1408) + 88 * (1:1408).^2), 1408);
    r_p2_interleaver = zeros(1, length(data));
    for i = (1:length(r_p2_interleaver))
        r_p2_interleaver(p2_interleaver(i) + 1) = data(i);
    end
end

function r_crc = crc(data)
    % data length is 1408
    % r_crc length is 1384
    % decode by 24A (24 last bits are crc)
    
    
    [~, err] = nrCRCDecode(data,"24A")
    % err_rate = sum(err)/(length(data) - 24);

    r_crc = data(1:end-24);
end

function r_turbo_and_interleaver = turbo_and_interleaver(data)
    % data length is 4236
    % r_turbo_and_interleaver length is 1408
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add noise
    % 
    % errors = randerr(1, length(data), 100);
    % data = xor(data, errors);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    S = data(1:1412);
    interlaced = data(1413:end);

    P1 = interlaced(1:2:end);
    P2 = interlaced(2:2:end);

    r_S = block_interleaver(S);
    r_P1 = block_interleaver(P1);
    r_P2 = block_interleaver(P2);
    
    % viterbi with P1 and then P2
    trellis = poly2trellis(4,[13 15], 13);

    code_1 = reshape([r_S; r_P1], 1, []);
    r_S_hat = vitdec(code_1, trellis, 20,'trunc','hard');    
    pi_S_hat = apply_p2_interleaver(r_S_hat);
    code_2 = reshape([pi_S_hat; r_P2], 1, []);
    r_turbo_and_interleaver = p2_interleaver(vitdec(code_2, trellis, 20,'trunc','hard'));
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % no error correction

    % r_turbo_and_interleaver = r_S;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % viterbi P1

    % code_1 = reshape([r_S; r_P1], 1, []);
    % r_turbo_and_interleaver = vitdec(code_1, trellis, 20,'trunc','hard');   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Turbo
    % n_iters = 100;
    % % S~~ = S
    % r_S_hat_hat = r_S;
    % metric_2 = [];
    % state_2 = [];
    % input_2 = [];
    % for i = (1:n_iters)
    %     % S~ = viterbi( S~~ , P1)
    %     code_1 = reshape([r_S_hat_hat; r_P1], 1, []);
    %     [r_S_hat, metric_1, state_1, input_1] = vitdec(code_1, trellis, 20, 'cont', 'unquant', metric_2, state_2, input_2);
    % 
    %     % S~~ = viterbi (pi(S~), P2)
    %     code_2 = reshape([apply_p2_interleaver(r_S_hat); r_P2], 1, []);
    %     [r_S_hat_hat, metric_2, state_2, input_2] = vitdec(code_2, trellis, 20,'cont','unquant', metric_1, state_1, input_1);
    %     r_S_hat_hat = p2_interleaver(r_S_hat_hat);
    % end
    % r_turbo_and_interleaver = r_S_hat_hat;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function r_block_interleaver = block_interleaver(data)
    % data length is 1412
    % r_block_interleaver legnth is 1408
    % matrix-interleaver (in by rows, out by columns)
    
    % add dummys
    dummys = [1 46 91 136 181 226 271 361 406 451 496 541 586 631 721 766 811 856 901 946 991 1081 1126 1171 1216 1261 1306 1351];
    for i = (1:length(dummys))
        value = dummys(i);
        data = [data(1:value - 1), NaN(1,1) ,data(value:end)];
    end

    % reshape to a matrix by columns
    matrix = reshape(data, 45, 32);

    % inverse the permutation
    % inverse of a cyclic permutation is the permutation flipped
    i_permutation = [1,17,9,25,5,21,13,29,3,19,11,27,7,23,15,31,2,18,10,26,6,22,14,30,4,20,12,28,8,24,16,32];
    new_matrix = zeros(45, 32);

    for i = (1:length(i_permutation))
        new_matrix(:,i) = matrix(:,i_permutation(i));
    end
    
    % reshape to a vector by rows
    r_block_interleaver = reshape(new_matrix.', 1, []);
     
    % remove dummys
    r_block_interleaver = r_block_interleaver(29:end);

    % remove 4 bits that turbo added
    r_block_interleaver = r_block_interleaver(1:end-4);
end

function r_cyclic_buffer = cyclic_buffer(data)
    % data size is 7200
    % removes cyclic buffer
    % r_cyclic_buffer is of length 4236

    r_cyclic_buffer = [data(4149:4236), data(1:4148)];
end

function r_xor_mask = xor_mask(data)
    % data size is 7200
    % xor with mask of length 7200
    mask = gen_gold_code(1600, 7200, 0x12345678).';
    r_xor_mask = bitxor(data, mask);
end

function seq = gen_gold_code(Nc, L, seed)
    % seq is column vector

    %%
    arguments
    Nc (1,1) = 1600
    L (1,1) = 7200
    seed (1,1) = 0x12345678
    end
    %%
    x1 = zeros(Nc + L + 31, 1);
    x2 = zeros(Nc + L + 31, 1);
    x2(1:32) = flip(fdec2bin(seed, 32));
    x1(1) = 1;
    for n = 1:(Nc + L)
    x1(n + 31) = xor(x1(n + 3), x1(n));
    x2(n + 31) = xor(xor(x2(n + 3), x2(n + 2)), xor(x2(n+1), x2(n)));
    end
    seq = xor(x1(Nc + 1:Nc+L), x2(Nc + 1:Nc+L));
end

function bits = fdec2bin(dec_num, n)
    bits = dec2bin(dec_num, n) - '0';
end




%% use later - TRYING TO HOP BETWEEN FREQUENCIES
function [fc_next] = hopping_between_freq(is_packet,fc_now, fc_list)
    % moving to next fc option
    idx = find(fc_list == fc_now);

    if ~is_packet
        if idx == length(fc_list)
            fc_next = fc_list(1);
        else
            fc_next = fc_list(idx + 1);
        end
    else
        fc_next = fc_now;
    end
end



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