clc
clear
close all

%% demod
fid = fopen('C:\Users\Student\Documents\final project\raw_bits_2.txt');
raw_hex = fread(fid,'*char').';

raw_bin = hexToBinaryVector(raw_hex, 4 * length(raw_hex));

r_xor_mask = xor_mask(raw_bin);
r_cyclic_buffer = cyclic_buffer(r_xor_mask);
r_turbo_and_intereaver = turbo_and_interleaver(r_cyclic_buffer);
r_crc = crc(r_turbo_and_intereaver);

proccesed_bin = r_crc(1:91*8);
proccesed_text = binaryVectorToHex(proccesed_bin);

%% check ber
fid = fopen('C:\Users\Student\Documents\final project\processed_bits12_Feb_2026_09_30_39_442.txt');
validation_hex = fread(fid,'*char').';
validation_bin = hexToBinaryVector(validation_hex, 4 * length(validation_hex));

ber = sum(validation_bin ~= proccesed_bin) / length(validation_bin)
[row, col] = find(validation_bin ~= proccesed_bin)

%% parse
% example = fopen('C:\Users\Student\Documents\final project\processed_bits12_Feb_2026_09_30_39_442.txt');
% example_hex = fread(example,'*char').';
% example_bin = hexToBinaryVector(example_hex);


lengths = [1 1 1 2 2 16 4 4 2 2 2 2 2 2 8 4 4 4 4 1 1 20 2];
offsets = [0 1 2 3 5 7 23 27 31 33 35 37 39 41 43 51 55 59 63 67 68 69 89];

for i = (1:length(offsets))
    current_bits = proccesed_bin(8*offsets(i) + 1:8*offsets(i) + 8*lengths(i));
    current_bits = current_bits.';
    if i == 1
        bit2int(current_bits, 8);
    end
    if i == 2
        % ...
    end
    if i == 3
        bit2int(current_bits, 8);
    end
    if i == 4
        bit2int(current_bits, 16);
    end
    if i == 5
        continue
        % ...
    end
    % PARTS ARE UINT8 AND ...
end

%% functions
function r_p2_interleaver = p2_interleaver(data)
    % data length is 1408
    % r_p2_interleaver length is 1408
    p2_interleaver = mod((43 * (1:1408) + 88 * (1:1408).^2), 1408);
    r_p2_interleaver = zeros(1, length(data));
    for i = (1:length(r_p2_interleaver))
        r_p2_interleaver(i) = data(p2_interleaver(i)+1);
    end
end

function r_crc = crc(data)
    % data length is 1408
    % r_crc length is 1384
    % decode by 24A (24 last bits are crc)
    
    
    [~, err] = nrCRCDecode(data,"24A");
    err_rate = sum(err)/(length(data) - 24);

    r_crc = data(1:end-24);
end

function r_turbo_and_interleaver = turbo_and_interleaver(data)
    % data length is 4236
    % r_turbo_and_interleaver length is 1408

    S = data(1:1412);
    interlaced = data(1413:end);

    P1 = interlaced(1:2:end);
    P2 = interlaced(2:2:end);

    r_S = block_interleaver(S);
    r_P1 = block_interleaver(P1);
    r_P2 = block_interleaver(P2);
    r_P2 = p2_interleaver(r_P2);

    % viterbi with P1
    trellis = poly2trellis(4,[13 15], 13);

    % interlacing
    code = reshape([r_S; r_P1], 1, []);

    r_turbo_and_interleaver = vitdec(code, trellis, 20,'trunc','hard');
    r_turbo_and_interleaver = r_S;
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

    r_cyclic_buffer = [data(4149:end), data(1:1184)];
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