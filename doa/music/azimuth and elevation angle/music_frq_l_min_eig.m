%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Multiple Signal Classification Algorithm for L-Shaped Array
%
% Description : MUltiple SIgnal Classification (MUSIC) Algorithm for L-Shaped Array
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 9 October 2016, Bandhit Suksiri,
%               Updated: 6 June    2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [norm_mat, sort_re_eig_val_vec] = music_frq_l_min_eig ( ...
	src_mat, n_sen_vec, d_sen, lamb, min_eig, src_frq, cen_frq, ang_az_rad_vec, ang_el_rad_vec)
	is_eulr_ang    = false;
    [n_tim, n_sen] = size(src_mat);
    
    cov_mat = zeros(n_sen, n_sen);
    for i = 1: 1: n_sen
        for j = 1: 1: n_sen
            cov_mat(i, j) = sum(src_mat(:, i) .* conj(src_mat(:, j)));
        end
    end
    cov_mat = (1 / n_tim) .* cov_mat;
    
    [eig_vec_mat, eig_val_mat] = eig(cov_mat);
    re_eig_val_vec             = zeros(n_sen, 1);
    for i = 1: 1: n_sen
        re_eig_val_vec(i, 1) = real(eig_val_mat(i, i));
    end
    [sort_re_eig_val_vec, idx_vec] = sort(re_eig_val_vec, 'descend');
    n_src = n_sen - sum(sort_re_eig_val_vec < min_eig);
    
    n_az_ang = size(ang_az_rad_vec, 1);
    n_el_ang = size(ang_el_rad_vec, 1);
    norm_mat = zeros(n_az_ang, n_el_ang);
    
    if (n_src > 0) && (n_src < n_sen)
        sort_eig_vec_mat = zeros(n_sen, n_sen);
        for i = 1: 1: n_sen
            sort_eig_vec_mat(:, i) = eig_vec_mat(:, idx_vec(i, 1));
        end
        e_mat     = sort_eig_vec_mat(:, n_src + 1: n_sen);
        sqr_e_mat = e_mat * e_mat';
        
        for i = 1: 1: n_az_ang
            for j = 1: 1: n_el_ang
                ang_rad_vec    = [ang_az_rad_vec(i, 1), ...
                                  ang_el_rad_vec(j, 1)];
                ster_vec       = l_fwd_ster_frq_fcn( ...
                    n_sen_vec, d_sen, ang_rad_vec, src_frq, cen_frq, lamb, is_eulr_ang);
                norm_mat(i, j) = 1 / (ster_vec' * sqr_e_mat * ster_vec);
            end
        end
        
        norm_mat = abs(norm_mat);
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%