%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : N. Tayem et al., L-Shape 2-Dimensional Arrival Angle Estimation with Propagator 
%               Method, IEEE, 2005.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 10 July 2017, Bandhit Suksiri,
%               Updated: 10 July 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_ang_x_rad_vec, out_ang_z_rad_vec] = doa_frq_n_src ( ...
    src_mat, n_sen, d_sen, lamb, n_src, src_frq, cen_frq)
    % n_sen_vec must be : [n_sen; 0; n_sen]!
    [n_tim, ~]  = size(src_mat);
    
    sub_x_mat   = src_mat(:, [1, n_sen + 1: 1: (2 * n_sen) - 2]).';
    sub_y_mat   = src_mat(:, [1, n_sen + 2: 1: (2 * n_sen) - 1]).';
    sub_z_mat   = src_mat(:, 1: 1: n_sen - 1).';
    sub_w_mat   = src_mat(:, 2: 1: n_sen).';
    
    r_xy_mat = zeros(2 * (n_sen - 1), 2 * (n_sen - 1));
    r_zw_mat = zeros(2 * (n_sen - 1), 2 * (n_sen - 1));
    for i = 1: 1: n_tim
        r_xy_mat = r_xy_mat + ([sub_x_mat(:, i); sub_y_mat(:, i)] * ...
                               [sub_x_mat(:, i); sub_y_mat(:, i)]');
        r_zw_mat = r_zw_mat + ([sub_z_mat(:, i); sub_w_mat(:, i)] * ...
                               [sub_z_mat(:, i); sub_w_mat(:, i)]');
    end
    r_xy_mat = r_xy_mat / n_tim;
    r_zw_mat = r_zw_mat / n_tim;
    
    e_xy_mat = r_xy_mat(:, 1: n_src);
    e_zw_mat = r_zw_mat(:, 1: n_src);
    j_xy_mat = r_xy_mat(:, n_src + 1: end);
    j_zw_mat = r_zw_mat(:, n_src + 1: end);
    
    p_csm_xy_mat = inv(e_xy_mat' * e_xy_mat) * e_xy_mat' * j_xy_mat; %#ok<MINV>
    p_csm_zw_mat = inv(e_zw_mat' * e_zw_mat) * e_zw_mat' * j_zw_mat; %#ok<MINV>
    
    p_csm_h_xy_mat = p_csm_xy_mat';
    p_csm_h_zw_mat = p_csm_zw_mat';
    
    p_1_xy_mat = p_csm_h_xy_mat(1: n_sen - 1 - n_src, :);
    p_1_zw_mat = p_csm_h_zw_mat(1: n_sen - 1 - n_src, :);
    p_3_xy_mat = p_csm_h_xy_mat(n_sen: end, :);
    p_3_zw_mat = p_csm_h_zw_mat(n_sen: end, :);
    
    [~, ohm_xy_mat] = eig(p_3_xy_mat * pinv(p_1_xy_mat));
    [~, ohm_zw_mat] = eig(p_3_zw_mat * pinv(p_1_zw_mat));
    
    eig_val_xy_vec  = diag(ohm_xy_mat);
    eig_val_zw_vec  = diag(ohm_zw_mat);
    [~, idx_xy_vec] = sort(abs(eig_val_xy_vec), 'descend');
    [~, idx_zw_vec] = sort(abs(eig_val_zw_vec), 'descend');
    eig_val_xy_vec  = eig_val_xy_vec(idx_xy_vec(1: n_src, 1), 1);
    eig_val_zw_vec  = eig_val_zw_vec(idx_zw_vec(1: n_src, 1), 1);
    
    out_ang_z_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                             angle(eig_val_xy_vec));
	out_ang_x_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                             angle(eig_val_zw_vec));
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%