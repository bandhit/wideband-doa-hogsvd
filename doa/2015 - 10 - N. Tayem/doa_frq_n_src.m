%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : N. Tayem et al., Two-Dimensional DOA Estimation Using Cross-Correlation Matrix With
%               L-Shaped Array, IEEE, 2015.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 8 July 2017, Bandhit Suksiri,
%               Updated: 8 July 2017, Bandhit Suksiri.
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
    
    r_xy_mat = zeros(n_sen - 1, n_sen - 1);
    r_zw_mat = zeros(n_sen - 1, n_sen - 1);
    for i = 1: 1: n_tim
        r_xy_mat = r_xy_mat + (sub_x_mat(:, i) * sub_y_mat(:, i)');
        r_zw_mat = r_zw_mat + (sub_z_mat(:, i) * sub_w_mat(:, i)');
    end
    r_xy_mat = r_xy_mat / n_tim;
    r_zw_mat = r_zw_mat / n_tim;
    
    e_xy_mat     = r_xy_mat + r_xy_mat';
    e_xy_3_2_mat = e_xy_mat((2 * n_src) + 1: n_sen - 1, n_src + 1: 2 * n_src);
    e_xy_2_1_mat = e_xy_mat(     n_src  + 1: 2 * n_src,         1:     n_src);
    e_xy_3_1_mat = e_xy_mat((2 * n_src) + 1: n_sen - 1,         1:     n_src);
    e_xy_1_2_mat = e_xy_mat(              1: n_src,     n_src + 1: 2 * n_src);
    e_zw_mat     = r_zw_mat + r_zw_mat';
    e_zw_3_2_mat = e_zw_mat((2 * n_src) + 1: n_sen - 1, n_src + 1: 2 * n_src);
    e_zw_2_1_mat = e_zw_mat(     n_src  + 1: 2 * n_src,         1:     n_src);
    e_zw_3_1_mat = e_zw_mat((2 * n_src) + 1: n_sen - 1,         1:     n_src);
    e_zw_1_2_mat = e_zw_mat(              1: n_src,     n_src + 1: 2 * n_src);
    
    phi_xy_1_mat = e_xy_3_2_mat * inv(e_xy_2_1_mat)';
    phi_xy_2_mat = e_xy_3_1_mat * inv(e_xy_1_2_mat)';
    phi_zw_1_mat = e_zw_3_2_mat * inv(e_zw_2_1_mat)';
    phi_zw_2_mat = e_zw_3_1_mat * inv(e_zw_1_2_mat)';
    
    phi_xy_mat = [phi_xy_1_mat; phi_xy_2_mat];
    phi_zw_mat = [phi_zw_1_mat; phi_zw_2_mat];
    
    o_xy_mat = phi_xy_mat(1: 2 * ((n_sen - 1) - (2 * n_src)) - 1, :);
    p_xy_mat = phi_xy_mat(2: 2 * ((n_sen - 1) - (2 * n_src)),     :);
    o_zw_mat = phi_zw_mat(1: 2 * ((n_sen - 1) - (2 * n_src)) - 1, :);
    p_zw_mat = phi_zw_mat(2: 2 * ((n_sen - 1) - (2 * n_src)),     :);
    
    [~, ohm_xy_mat] = eig(pinv(o_xy_mat) * p_xy_mat);
    [~, ohm_zw_mat] = eig(pinv(o_zw_mat) * p_zw_mat);
    
    out_ang_z_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                             angle(diag(ohm_xy_mat)));
	out_ang_x_rad_vec = acos(((lamb * cen_frq) / (2 * pi * d_sen * src_frq)) * ...
                             angle(diag(ohm_zw_mat)));
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%