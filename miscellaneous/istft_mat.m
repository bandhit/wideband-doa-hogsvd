%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Inverse Short-time Fourier Transform Function 
%
% Description : Inverse Short-time Fourier Transform (ISTFT) with Multiple Vector Representation
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 12 November 2016, Bandhit Suksiri,
%               Updated:  6 June     2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [syn_mat, tim_vec] = istft_mat (stft_mat3, sam_frq, win_fcn, n_fft, hop_size)
    [~, n_tim_fft, n_sam] = size(stft_mat3);
    n_tim                 = n_fft + (n_tim_fft - 1) * hop_size;
    win_vec               = feval(win_fcn, n_fft);
    tim_vec               = (0: 1: n_tim - 1)' ./ sam_frq;
    
    syn_mat = zeros(n_tim, n_sam);
    for j = 1: 1: n_sam
        if rem(n_fft, 2)
            for i = 0: hop_size: (hop_size * (n_tim_fft - 1))
                fft_vec = stft_mat3(:, 1 + i / hop_size, j);
                fft_mat = [fft_vec; conj(fft_vec(end: -1: 2))];
                syn_vec = real(ifft(fft_mat));
                syn_vec = win_vec .* syn_vec;
                syn_mat(i + 1: i + n_fft, j) = syn_mat(i + 1: i + n_fft, j) + syn_vec;
            end
        else
            for i = 0: hop_size: (hop_size * (n_tim_fft - 1))
                fft_vec = stft_mat3(:, 1 + i / hop_size, j);
                fft_mat = [fft_vec; conj(fft_vec(end - 1: -1: 2))];
                syn_vec = real(ifft(fft_mat));
                syn_vec = win_vec .* syn_vec;
                syn_mat(i + 1: i + n_fft, j) = syn_mat(i + 1: i + n_fft, j) + syn_vec;
            end
        end
    end
    
    k       = sum(win_vec .^ 2);
    syn_mat = (syn_mat .* hop_size) ./ k;
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%