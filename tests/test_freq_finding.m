

nfft = 1024; % Number of FFT points

% ---------create zadof chu------- %
Nzc = 601;

root = 600;
seq_4 = zadof_ofdm(Nzc,root);

root = 147;
seq_6 = zadof_ofdm(Nzc,root);

t = 0:1/fs:length(seq_4)/fs-1/fs;

freq_offset = 1001.21;

seq_4_shift = seq_4 .* exp(-1j*2*pi*freq_offset.*t);


freq_offset = kay(seq_4,seq_4_shift,fs)


function freq_offset = kay(zadof,pilot,fs)
    L = length(pilot);
    k = (1:L-1);
    z = conj(zadof) .* pilot;

    z2 = z(1:end - 1);
    z1 = z(2:end);
    
    w = (3/2) * (L / (L^2 - 1)) * (1 - ((2 * k - L) / L).^2);

    freq_offset = sum(w .* (angle(z1 .* conj(z2)))) * fs / (2 * pi);
end


function seq_freq = zadof_ofdm(Nzc,root)
    
    m = 0:(Nzc-1);
    seq = exp((-1j*pi*root*m.*(m+1)/Nzc));

    seq = [zeros(1,212),seq,zeros(1,(1024-813))];
  
    seq_freq = ifft(ifftshift(seq));

end


