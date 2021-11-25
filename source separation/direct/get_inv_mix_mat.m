%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Inverse mixture matrix of mixed-signal for L-Shaped Array
%
% Description : Inverse mixture matrix of mixed-signal for L-Shaped Array,
%               proof-of-concept
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 14 January 2019, Bandhit Suksiri,
%               Updated: 15 January 2019, Bandhit Suksiri.
%
% Copyright 2019,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inv_mix_mat3 = get_inv_mix_mat (ang_rad_mat, frq_vec, n_sen, d_sen, lamb, cen_frq)
    is_eulr_ang = true;
    n_frq       = size(frq_vec, 1);
    n_src       = size(ang_rad_mat, 1);
    n_sen_vec   = [n_sen; 0; n_sen - 1];
    
    inv_mix_mat3 = zeros(n_src, 2 * n_sen, n_frq);
    for i_frq = 1: 1: n_frq
        tmp_src_frq   = frq_vec(i_frq, 1);
        calc_a_fi_mat = l_fwd_ster_frq_fcn( ...
            n_sen_vec, d_sen, ang_rad_mat, tmp_src_frq, cen_frq, lamb, is_eulr_ang);
        
%         inv_mix_mat3(:, :, i_frq) = pinv(calc_a_fi_mat);
        
        calc_a_fi_x_mat = calc_a_fi_mat(2: 1: n_sen + 1, :);
        calc_a_fi_z_mat = calc_a_fi_mat([1, n_sen + 2: 1: 2 * n_sen], :);
        
        inv_mix_mat3(:, 2: 1: n_sen + 1,              i_frq) = pinv(calc_a_fi_x_mat);
        inv_mix_mat3(:, [1, n_sen + 2: 1: 2 * n_sen], i_frq) = pinv(calc_a_fi_z_mat);
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
