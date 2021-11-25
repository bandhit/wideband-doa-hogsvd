%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for Uniform Linear Array
%
% Description : J. Cao et al., DOA Estimation for Wideband Sources Using Cross Correlation 
%               Transformation, IEEE, 2010.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 19 June 2017, Bandhit Suksiri,
%               Updated: 21 June 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [norm_vec, sort_re_eig_val_vec] = doa_ula_n_src ( ...
    src_mat3, d_sen, lamb, n_src, ...
    idx_ref, ang_rad_vec, axis_sel)
    [n_frq, n_tim, n_sen] = size(src_mat3);
    
    r_fi_f0_mat = zeros(n_sen, n_sen, n_frq);
    for i = 1: 1: n_frq
        for j = 1: 1: n_tim
            r_fi_f0_mat(:, :, i) = r_fi_f0_mat(:, :, i) ...
                                 + (squeeze(src_mat3(i, j, :)) * squeeze(src_mat3(idx_ref, j, :))');
        end
    end
    r_fi_f0_mat = r_fi_f0_mat / n_tim;
    
    p_mat = zeros(n_sen, n_sen);
    for j = 1: 1: n_tim
        p_mat(:, :) = p_mat(:, :) ...
                    + (squeeze(src_mat3(idx_ref, j, :)) * squeeze(src_mat3(idx_ref, j, :))');
    end
    p_mat = p_mat / n_tim;
    
    r_mat = zeros(n_sen, n_sen);
    for i = 1: 1: n_frq
        [u_mat, ~, v_mat] = svd(r_fi_f0_mat(:, :, i) * p_mat');
        t_mat             = v_mat * u_mat';
        
        r_fi_fi_mat = zeros(n_sen, n_sen);
        for j = 1: 1: n_tim
            r_fi_fi_mat(:, :) = r_fi_fi_mat(:, :) ...
                              + (squeeze(src_mat3(i, j, :)) * squeeze(src_mat3(i, j, :))');
        end
        r_fi_fi_mat = r_fi_fi_mat / n_tim;
        r_mat = r_mat + (t_mat * r_fi_fi_mat * t_mat');
    end
    
    [eig_vec_mat, eig_val_mat] = eig(r_mat);
    re_eig_val_vec             = zeros(n_sen, 1);
    for i = 1: 1: n_sen
        re_eig_val_vec(i, 1) = real(eig_val_mat(i, i));
    end
    [sort_re_eig_val_vec, idx_vec] = sort(re_eig_val_vec, 'descend');
    
    n_ang = size(ang_rad_vec, 1);
    if axis_sel == 'x'
        ang_x_rad_vec = ang_rad_vec;
        ang_y_rad_vec = ndef_ang * ones(n_ang, 1);
        ang_z_rad_vec = ndef_ang * ones(n_ang, 1);
    elseif axis_sel == 'y'
        ang_x_rad_vec = ndef_ang * ones(n_ang, 1);
        ang_y_rad_vec = ang_rad_vec;
        ang_z_rad_vec = ndef_ang * ones(n_ang, 1);
    elseif axis_sel == 'z'
        ang_x_rad_vec = ndef_ang * ones(n_ang, 1);
        ang_y_rad_vec = ndef_ang * ones(n_ang, 1);
        ang_z_rad_vec = ang_rad_vec;
    end
    norm_vec = zeros(n_ang, 1);
    
    if (n_src > 0) && (n_src < n_sen)
        sort_eig_vec_mat = zeros(n_sen, n_sen);
        for i = 1: 1: n_sen
            sort_eig_vec_mat(:, i) = eig_vec_mat(:, idx_vec(i, 1));
        end
        e_mat     = sort_eig_vec_mat(:, n_src + 1: n_sen);
        sqr_e_mat = e_mat * e_mat';
        
        for i = 1: 1: n_ang
            tmp_ang_rad_vec = [ang_x_rad_vec(i, 1), ...
                               ang_y_rad_vec(i, 1), ...
                               ang_z_rad_vec(i, 1)];
            ster_vec        = ula_fwd_ster_fcn(n_sen, d_sen, tmp_ang_rad_vec, lamb, axis_sel);
            norm_vec(i, 1)  = 1 / (ster_vec' * sqr_e_mat * ster_vec);
        end
        
        norm_vec = abs(norm_vec);
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%