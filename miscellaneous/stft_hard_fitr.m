%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Short-time Fourier Transform Post-function, Frequency Selection
%
% Description : Short-time Fourier Transform Post-function, Frequency Selection
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

function [post_stft_mat3, post_frq_vec, head_idx, tail_idx] = stft_hard_fitr ( ...
    stft_mat3, frq_vec, min_frq_vec, max_frq_vec)
    idx_vec = find((frq_vec >= min_frq_vec) & (frq_vec <= max_frq_vec));
    if ~isempty(idx_vec)
        head_idx       = idx_vec(1, 1);
        tail_idx       = idx_vec(end, 1);
        post_frq_vec   = frq_vec(head_idx: tail_idx, 1);
        post_stft_mat3 = stft_mat3(head_idx: tail_idx, :, :);
    else
        head_idx       = 1;
        tail_idx       = size(frq_vec, 1);
        post_frq_vec   = frq_vec;
        post_stft_mat3 = stft_mat3;
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
