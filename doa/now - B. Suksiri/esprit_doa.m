%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : ESPRIT DOA Estimation
%
% Description : ESPRIT DOA Estimation, IEEE
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

function [calc_ang_x_rad_vec, eig_vec_x_mat, eig_val_x_vec] = esprit_doa ...
        (v_x_mat, fc_src_frq, cen_frq, lamb, d_sen, n_rank_fc, is_pinv, is_eig)
    
    v_x_1_mat = v_x_mat(1: end - 1, end - n_rank_fc + 1: end);
    v_x_2_mat = v_x_mat(2: end,     end - n_rank_fc + 1: end);
    
    if is_pinv == true
        ang_x_mat = pinv(v_x_1_mat) * v_x_2_mat;
    else
        ang_x_mat = lsqminnorm(v_x_1_mat, v_x_2_mat);
    end
    
    if is_eig == true
        [eig_vec_x_mat, eig_val_x_mat] = eig(ang_x_mat);
    else
        [eig_vec_x_mat, eig_val_x_mat] = schur(ang_x_mat, 'complex');
    end
    eig_val_x_vec = diag(eig_val_x_mat);
    
    calc_ang_x_rad_vec = esprit_ang_fcn(eig_val_x_vec, fc_src_frq, cen_frq, lamb, d_sen);
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%