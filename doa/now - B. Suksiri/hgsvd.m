%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Higher-Order Generalized Singular Value Decomposition (hgsvd)
%
% Description : https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0028072
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

function [v_mat, u_mat3, s_mat3, l_mat] = hgsvd (full_a_mat, n_rank_gsvd, is_qr_schr, is_pinv, is_get_u)
    
    if is_pinv == true
        inv_fcn = @pinv;
    else
        inv_fcn = @inv;
    end
    
    [n_row, n_col] = size(full_a_mat);
    n_sam          = n_row / n_col;
    
    
    if is_qr_schr == true
        [full_q_mat, r_vec] = qr(full_a_mat, 0);
    else
        full_q_mat = full_a_mat;
    end
    
    t_mat = zeros(n_col, n_col);
    for i_sam = 1: 1: n_sam
        q_mat = full_q_mat(((i_sam - 1) * n_col) + 1: i_sam * n_col, :);
        t_mat = t_mat + inv_fcn(q_mat' * q_mat);
    end
    t_mat = t_mat ./ n_sam;
    
    if is_qr_schr == true
        [z_mat, l_mat] = schur(t_mat, 'real');
    else
        [z_mat, l_mat] = eig(t_mat);
    end
    
    l_mat        = real(l_mat);
    [~, idx_vec] = sort(diag(l_mat), 'descend');

    l_mat = l_mat(idx_vec, idx_vec);
    z_mat = z_mat(:, idx_vec);
    
    l_mat = l_mat(end - n_rank_gsvd + 1: end, end - n_rank_gsvd + 1: end);
    z_mat = z_mat(:, end - n_rank_gsvd + 1: end);
    
    if is_qr_schr == true
        v_mat = r_vec' * z_mat;
    else
        v_mat = z_mat;
    end
    
    s_mat3 = zeros(n_rank_gsvd, n_rank_gsvd, n_sam);
    u_mat3 = zeros(n_col,       n_rank_gsvd, n_sam);
    if is_get_u == true
        for i_sam = 1: 1: n_sam
            q_mat  = full_q_mat(((i_sam - 1) * n_col) + 1: i_sam * n_col, :);
            qz_mat = q_mat * z_mat;

            s_mat = diag(sqrt(sum(abs(qz_mat) .^ 2,  1)));
            u_mat = qz_mat / s_mat;

            s_mat3(:, :, i_sam) = s_mat;
            u_mat3(:, :, i_sam) = u_mat;
            
%             u_mat(:, 4: end)' * u_mat(:, 4: end)
%             u_mat(:, 1: 3)' * u_mat(:, 1: 3)
%             u_mat(:, 4: end)' * u_mat(:, 1: 3)
%             u_mat(:, 1: 3)' * u_mat(:, 4: end)
%             u_mat(:, 4: end)' * u_mat(:, 1: end)
%             u_mat(:, 1: 3)' * u_mat(:, 1: end)
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%