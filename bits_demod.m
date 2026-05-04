
raw_bin = ??;

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

%% FUNCTIONS
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
    
    
    [~, err] = nrCRCDecode(data,"24A");
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