%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Get First N peaks of Angle Function on 2 Dimensional Information
%
% Description : Get first maximal N peaks
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

function [ang_rad_x_vec, ang_rad_z_vec] = get_peak_2d (abs_norm_mat, ang_rad_vec, n_src)
    n_ang = size(ang_rad_vec, 1);
    
    x_pks_vec    = [];
    x_locs_x_vec = [];
    x_locs_z_vec = [];
    for i_ang = 1: n_ang
        [cnt_x_pks_cvec, cnt_x_locs_z_cvec] = findpeaks(abs_norm_mat(i_ang, :));
        
        x_pks_vec    = [x_pks_vec;    cnt_x_pks_cvec.']; %#ok<AGROW>
        x_locs_x_vec = [x_locs_x_vec; i_ang .* ones(size(cnt_x_pks_cvec, 2), 1)]; %#ok<AGROW>
        x_locs_z_vec = [x_locs_z_vec; cnt_x_locs_z_cvec.']; %#ok<AGROW>
    end
    
    [~, x_idx_vec] = sort(x_pks_vec, 'descend');
    x_locs_x_vec   = x_locs_x_vec(x_idx_vec, 1);
    x_locs_z_vec   = x_locs_z_vec(x_idx_vec, 1);
    
    n_x_peak   = size(x_pks_vec, 1);
    locs_x_vec = [];
    locs_z_vec = [];
    i_src      = 1;
    for i_x_peak = 1: n_x_peak
        x_locs_x = x_locs_x_vec(i_x_peak, 1);
        x_locs_z = x_locs_z_vec(i_x_peak, 1);
        
        if (x_locs_x > 1) && (x_locs_x < n_ang) && (x_locs_z > 1) && (x_locs_z < n_ang)
            raw_11 = abs_norm_mat(x_locs_x - 1, x_locs_z - 1);
            raw_21 = abs_norm_mat(x_locs_x,     x_locs_z - 1);
            raw_31 = abs_norm_mat(x_locs_x + 1, x_locs_z - 1);
            raw_12 = abs_norm_mat(x_locs_x - 1, x_locs_z);
            raw_22 = abs_norm_mat(x_locs_x,     x_locs_z);
            raw_32 = abs_norm_mat(x_locs_x + 1, x_locs_z);
            raw_13 = abs_norm_mat(x_locs_x - 1, x_locs_z + 1);
            raw_23 = abs_norm_mat(x_locs_x,     x_locs_z + 1);
            raw_33 = abs_norm_mat(x_locs_x + 1, x_locs_z + 1);
            if (raw_22 > raw_11) && ...
               (raw_22 > raw_21) && ...
               (raw_22 > raw_31) && ...
               (raw_22 > raw_12) && ...
               (raw_22 > raw_32) && ...
               (raw_22 > raw_13) && ...
               (raw_22 > raw_23) && ...
               (raw_22 > raw_33)
                locs_x_vec = [locs_x_vec; x_locs_x]; %#ok<AGROW>
                locs_z_vec = [locs_z_vec; x_locs_z]; %#ok<AGROW>
                
                i_src = i_src + 1;
                
                if i_src > n_src
                    break;
                end
            end
        end
    end
    
    post_n_src = size(locs_x_vec, 1);
    if post_n_src > 0
        if post_n_src > n_src
            post_n_src = n_src;
        end
        ang_rad_x_vec = ang_rad_vec(locs_x_vec(1: post_n_src, 1), 1);
        ang_rad_z_vec = ang_rad_vec(locs_z_vec(1: post_n_src, 1), 1);

        n_lost = n_src - post_n_src;
        if n_lost > 0
            ang_rad_x_vec = [zeros(n_lost, 1); ang_rad_x_vec];
            ang_rad_z_vec = [zeros(n_lost, 1); ang_rad_z_vec];
        end
    else
        ang_rad_x_vec = zeros(n_src, 1);
        ang_rad_z_vec = zeros(n_src, 1);
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%