%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Synthesis Sources Function
%
% Description : Synthesis Sources 
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 16 November 2016, Bandhit Suksiri,
%               Updated:  6 June     2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [syn_mat, tim_istft_vec, syn_sfft_mat3] = l_syn_fcn ( ...
        ang_az_rad_vec, ang_el_rad_vec, n_sen_vec, d_sen, cen_frq, lamb, ...
        sfft_mat3, frq_vec, sam_frq, win_fcn, n_fft, hop_size, min_frq, max_frq)
    is_eulr_ang = false;
    ang_rad_vec = [ang_az_rad_vec(:, 1), ...
                   ang_el_rad_vec(:, 1)];
	
    [~, ~, head_idx, tail_idx] = stft_hard_fitr(sfft_mat3, frq_vec, min_frq, max_frq);
    [n_frq_fft, n_tim_fft, ~]  = size(sfft_mat3);
    
    syn_sfft_mat3 = zeros(n_frq_fft, n_tim_fft, size(ang_rad_vec, 1));
    for i = head_idx: 1: tail_idx
        src_fft_mat = squeeze(sfft_mat3(i, :, :))';
        src_frq     = frq_vec(i, 1);
        ster_mat    = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, ang_rad_vec, ...
                                         src_frq, cen_frq, lamb, ...
                                         is_eulr_ang) .^ -1;
        syn_sfft_mat3(i, :, :) = (pinv(ster_mat) * src_fft_mat)';
    end
    [syn_mat, tim_istft_vec] = istft_mat(syn_sfft_mat3, sam_frq, win_fcn, n_fft, hop_size);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%