%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : MUSIC DOA Estimation
%
% Description : MUSIC DOA Estimation, IEEE
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 14 March 2019, Bandhit Suksiri,
%               Updated: 14 March 2019, Bandhit Suksiri.
%
% Copyright 2019,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function norm_mat = music_doa_sig (v_mat, n_sen_vec, fc_src_frq, cen_frq, lamb, d_sen, n_rank_fc, ang_x_rad_vec, ang_z_rad_vec, ster_idx_cvec)
    
%     sqr_e_x_mat = v_mat(:, end - n_rank_fc: -1: 1) * v_mat(:, end - n_rank_fc: -1: 1)';
    sqr_s_x_mat = v_mat(:, end - n_rank_fc + 1: end) * v_mat(:, end - n_rank_fc + 1: end)';
    n_sen = size(v_mat, 1);
    
    n_ang_x  = size(ang_x_rad_vec, 1);
    n_ang_z  = size(ang_z_rad_vec, 1);
    norm_mat = zeros(n_ang_x, n_ang_z);
    for i_ang_x = 1: 1: n_ang_x
        for i_ang_z = 1: 1: n_ang_z
            iter_ang_rad_vec = [ang_x_rad_vec(i_ang_x, 1), ndef_ang, ang_z_rad_vec(i_ang_z, 1)];

            ster_vec = l_fwd_ster_frq_fcn(n_sen_vec, d_sen, iter_ang_rad_vec, fc_src_frq, cen_frq, lamb, true);
            ster_vec = ster_vec(ster_idx_cvec, :);
            
            norm_mat(i_ang_x, i_ang_z) = (ster_vec' * ster_vec) ./ ((ster_vec' * ster_vec) - (ster_vec' * sqr_s_x_mat * ster_vec));
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%