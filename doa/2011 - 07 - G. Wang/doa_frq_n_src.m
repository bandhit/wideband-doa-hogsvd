%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Direction of Arrival (DOA) Estimation for L-Shaped Array
%
% Description : G. Wang et al., Computationally Efficient Subspace-Based Method for Two-Dimensional 
%               Direction Estimation with L-Shaped Array, IEEE, 2011.
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 8 June 2017, Bandhit Suksiri,
%               Updated: 8 June 2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function norm_mat = doa_frq_n_src ( ...
    src_mat, n_sen, d_sen, lamb, n_src, src_frq, cen_frq, ...
    ang_x_rad_vec, ang_z_rad_vec)
    % n_sen_vec must be : [n_sen; 0; n_sen - 1]!
	is_eulr_ang = true;
    [n_tim, ~]  = size(src_mat);
    n_sen_vec   = [n_sen; 0; n_sen - 1];
    x_mat       = src_mat(:, 2: 1: n_sen + 1).';
    z_mat       = src_mat(:, [1, n_sen + 2: 1: 2 * n_sen]).';
    y_mat       = [x_mat; z_mat];
    
    r_yy_mat = zeros(2 * n_sen, 2 * n_sen);
    for i = 1: 1: n_tim
        r_yy_mat = r_yy_mat + (y_mat(:, i) * y_mat(:, i)');
    end
    r_yy_mat = r_yy_mat / n_tim;
    
    g1_mat = r_yy_mat(        1:     n_src, :);
    g2_mat = r_yy_mat(n_src + 1: 2 * n_sen, :);
    p_mat  = inv(g1_mat * g1_mat') * (g1_mat * g2_mat'); %#ok<MINV>
    q_mat  = [p_mat.', - eye((2 * n_sen) - n_src)].';
    pi_mat = q_mat ...
           * (eye((2 * n_sen) - n_src) - ...
              (p_mat' * ...
               inv((p_mat * p_mat') + eye(n_src)) * ...
               p_mat)) ...
           * q_mat'; %#ok<MINV>
    
    n_ang_x  = size(ang_x_rad_vec, 1);
    n_ang_z  = size(ang_z_rad_vec, 1);
    norm_mat = zeros(n_ang_x, n_ang_z);
    
    for i = 1: 1: n_ang_x
        for j = 1: 1: n_ang_z
            ang_rad_vec    = [ang_x_rad_vec(i, 1), ...
                              ndef_ang, ...
                              ang_z_rad_vec(j, 1)];
            ster_vec       = l_fwd_ster_frq_fcn( ...
                n_sen_vec, d_sen, ang_rad_vec, src_frq, cen_frq, lamb, is_eulr_ang);
            norm_mat(i, j) = ster_vec' * pi_mat * ster_vec;
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%