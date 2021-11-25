%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Short-time Fourier Transform Function 
%
% Description : Short-time Fourier Transform (STFT) with Multiple Vector Representation (Matrix)
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 12 November 2016, Bandhit Suksiri,
%               Updated: 18 May      2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [stft_mat3, frq_vec, tim_vec] = stft_mat (src_mat, sam_frq, win_fcn, win_size, n_fft, ...
                                                   hop_size)
    [n_tim, n_sam] = size(src_mat);
    win_vec        = feval(win_fcn, win_size);
    n_frq_fft      = ceil((n_fft + 1) / 2);
    n_tim_fft      = fix((n_tim - win_size) / hop_size) + 1;
    frq_vec        = (sam_frq / n_fft) .* (0: 1: n_frq_fft - 1)';
    tim_vec        = (win_size / 2: ...
                      hop_size: ...
                      win_size / 2 + (n_tim_fft - 1) * hop_size)' ...
                   / sam_frq;
    
    head_hop_vec = zeros(n_tim_fft, 1);
    tail_hop_vec = zeros(n_tim_fft, 1);
    idx_hop      = 1;
    for i = 1: 1: n_tim_fft
        head_hop_vec(i, 1) = idx_hop;
        tail_hop_vec(i, 1) = idx_hop + win_size - 1;
        idx_hop            = idx_hop + hop_size;
    end
    
    stft_mat3 = zeros(n_frq_fft, n_tim_fft, n_sam);
    for k = 1: 1: n_sam
        src_vec = src_mat(:, k);
        for j = 1: 1: n_tim_fft
            conv_vec           = src_vec(head_hop_vec(j, 1): tail_hop_vec(j, 1), 1) .* win_vec;
            fft_vec            = fft(conv_vec, n_fft);
            stft_mat3(:, j, k) = fft_vec(1: n_frq_fft, 1);
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%