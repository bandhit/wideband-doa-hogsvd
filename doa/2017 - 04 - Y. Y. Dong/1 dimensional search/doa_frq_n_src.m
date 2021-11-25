%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : Y. Y. Dong et al., 2-D DOA Estimation for L-Shaped Array With Array Aperture and 
%               Snapshots Extension Techniques, IEEE, 2017.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 5 June 2017, Bandhit Suksiri,
%               Updated: 8 June 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out_ang_x_rad_vec, out_ang_z_rad_vec, sort_re_r_eig_val_vec] = ...
    doa_frq_n_src (src_mat, n_sen, d_sen, lamb, n_src, src_frq, cen_frq, ...
    ang_x_rad_vec)
    % n_sen_vec must be : [n_sen; 0; n_sen - 1]!
	is_eulr_ang = true;
    [n_tim, ~]  = size(src_mat);
    x_mat       = src_mat(:, 2: 1: n_sen + 1).';
    z_mat       = src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    x_bar_mat   = src_mat(:, n_sen + 1: -1: 2).';
    
    r_xz_mat     = zeros(n_sen, n_sen);
    r_bar_xz_mat = zeros(n_sen, n_sen);
    for i = 1: 1: n_tim
        r_xz_mat     = r_xz_mat     + (x_mat(:, i)     * z_mat(:, i)');
        r_bar_xz_mat = r_bar_xz_mat + (x_bar_mat(:, i) * z_mat(:, i).');
    end
    r_xz_mat     = r_xz_mat     / n_tim;
    r_bar_xz_mat = r_bar_xz_mat / n_tim;
    
    y1_mat     = r_xz_mat(:, 1: n_sen - 1);
    y2_mat     = r_xz_mat(:, 2: n_sen);
    y1_bar_mat = r_bar_xz_mat(:, 1: n_sen - 1);
    y2_bar_mat = r_bar_xz_mat(:, 2: n_sen);
    
    y_mat  = [y1_mat, y2_bar_mat;
              y2_mat, y1_bar_mat];
    
	r_yy_mat = y_mat * y_mat';
    
    [r_eig_vec_mat, r_eig_val_mat] = eig(r_yy_mat);
    re_r_eig_val_vec               = zeros(2 * n_sen, 1);
    for i = 1: 1: 2 * n_sen
        re_r_eig_val_vec(i, 1) = real(r_eig_val_mat(i, i));
    end
    [sort_re_r_eig_val_vec, sort_re_r_eig_idx_vec] = sort(re_r_eig_val_vec, 'descend');
    
    
    if (n_src > 0) && (n_src < 2 * n_sen)
        sort_eig_vec_mat = zeros(2 * n_sen, 2 * n_sen);
        for i = 1: 1: 2 * n_sen
            sort_eig_vec_mat(:, i) = r_eig_vec_mat(:, sort_re_r_eig_idx_vec(i, 1));
        end
        u_w_mat = sort_eig_vec_mat(:, n_src + 1: 2 * n_sen);
        
        a_x_fcn = @(ang_x_rad) fwd_ster_frq_fcn(d_sen * [(1: 1: n_sen)', ...
                                                         zeros(n_sen, 1), ...
                                                         zeros(n_sen, 1)], ...
                                                [ang_x_rad, ndef_ang, ndef_ang], ...
                                                src_frq, cen_frq, lamb, is_eulr_ang);
        f_fcn   = @(ang_x_rad) (kron(eye(2), a_x_fcn(ang_x_rad))' * ...
                                (u_w_mat * u_w_mat') * ...
                                kron(eye(2), a_x_fcn(ang_x_rad)));
        n_f     = 2; % always, from above
        
        n_ang_x   = size(ang_x_rad_vec, 1);
        f_out_vec = zeros(n_ang_x, 1);
        for i = 1: 1: n_ang_x
            f_out_vec(i, 1) = 1 / det(f_fcn(ang_x_rad_vec(i, 1)));
        end
        f_out_vec             = real(f_out_vec);
        [pks_vec, locs_vec]   = findpeaks(f_out_vec);
        [~, sort_pks_idx_vec] = sort(pks_vec, 'descend');
        n_pks                 = size(sort_pks_idx_vec, 1);
        if (n_pks < n_src)
            n_src = n_pks;
        end
        out_ang_x_rad_vec = ang_x_rad_vec(locs_vec(sort_pks_idx_vec(1 : n_src, 1), 1), 1);
        out_ang_z_rad_vec = zeros(n_src, 1);
        
        for i = 1: 1: n_src
            f_mat                          = f_fcn(out_ang_x_rad_vec(i, 1));
            [f_eig_vec_mat, f_eig_val_mat] = eig(f_mat);
            f_eig_val_vec                  = zeros(n_f, 1);
            for j = 1: 1: n_f
                f_eig_val_vec(j, 1) = real(f_eig_val_mat(j, j));
            end
            [~, sort_f_eig_idx_vec] = sort(f_eig_val_vec, 'ascend');
            sort_f_eig_vec_mat = f_eig_vec_mat(:, sort_f_eig_idx_vec(:, 1));
            out_ang_z_rad_vec(i, 1) = acos(angle(+ sort_f_eig_vec_mat(1, 1) / ...
                                                   sort_f_eig_vec_mat(2, 1)) * ...
                                           ((lamb * cen_frq) / (2 * pi * d_sen * src_frq)));
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%