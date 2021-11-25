%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Expectation Maximization Algorithm for Gaussian Mixture Model
%
% Description : Expectation Maximization Algorithm for Gaussian Mixture Model (GMM-EM)
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 22 November 2016, Bandhit Suksiri,
%               Updated: 26 May      2017, Bandhit Suksiri.
%
% Copyright 2016 - 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_mat, llh_vec, w_vec, mu_mat, sigma_mat3] = gmm_em (x_mat, n_k, max_iter, reg_val, tol)
    r_mat         = init_r_kmeans(x_mat, n_k);
    llh_vec       = zeros(max_iter, 1);
    llh_vec(1, 1) = -inf;
    
    for i = 2: 1: max_iter
        [w_vec, mu_mat, sigma_mat3] = m_step(x_mat, r_mat, reg_val);
        [r_mat, llh]                = e_step(x_mat, w_vec, mu_mat, sigma_mat3);

        llh_vec(i, 1) = llh;
        if abs(llh - llh_vec(i - 1, 1)) < (tol * abs(llh))
            llh_vec = llh_vec(1: i);
            break;
        end
    end
end

function r_mat = init_r_kmeans (x_mat, n_k)
    [n_dim, n_sam] = size(x_mat);
    if (n_dim >= 1) && (n_dim <= 2)
        r_mat   = zeros(n_sam, n_k);
        idx_vec = kmeans(x_mat', n_k, 'MaxIter', 1e3, 'Replicates', 1e1);
        for i = 1: 1: n_sam
            r_mat(i, idx_vec(i, 1)) = 1;
        end
    else
        r_mat = init_r_rand(x_mat, n_k);
    end
end

function r_mat = init_r_rand (x_mat, n_k)
    n_sam   = size(x_mat, 2);
    r_mat   = zeros(n_sam, n_k);
    idx_vec = randi([1, n_k], n_sam, 1);
    for i = 1: 1: n_sam
        r_mat(i, idx_vec(i, 1)) = 1;
    end
end

function [w_vec, mu_mat, sigma_mat3] = m_step (x_mat, r_mat, reg_val)
    [n_dim, n_sam] = size(x_mat);
    n_k            = size(r_mat, 2);
    
    m_col = sum(r_mat, 1);
    m_vec = m_col';
    w_vec = m_vec ./ n_sam;
    
    mu_mat = zeros(n_dim, n_k);
    for i = 1: 1: n_k
        for j = 1: 1: n_sam
            mu_mat(:, i) = mu_mat(:, i) ...
                         + (r_mat(j, i) .* x_mat(:, j));
        end
        mu_mat(:, i) = mu_mat(:, i) ./ m_vec(i, 1);
    end
    
    reg_mat    = reg_val * eye(n_dim);
    sigma_mat3 = zeros(n_dim, n_dim, n_k);
    for i = 1: 1: n_k
        for j = 1: 1: n_sam
            x_mu_vec            = x_mat(:, j) - mu_mat(:, i);
            sigma_mat3(:, :, i) = sigma_mat3(:, :, i) ...
                                + (r_mat(j, i) .* (x_mu_vec * x_mu_vec'));
        end
        sigma_mat3(:, :, i) = (sigma_mat3(:, :, i) ./ m_vec(i, 1)) + reg_mat;
    end
end

function [new_r_mat, llh] = e_step (x_mat, w_vec, mu_mat, sigma_mat3)
    n_sam     = size(x_mat, 2);
    n_k       = size(w_vec, 1);
    
    new_r_mat = zeros(n_sam, n_k);
    for i = 1: 1: n_k
        new_r_mat(:, i) = w_vec(i, 1) * mgm(x_mat, mu_mat(:, i), sigma_mat3(:, :, i));
    end
    
    sum_gmm_vec = sum(new_r_mat, 2);
    llh         = sum(log(sum_gmm_vec)) ./ n_sam;
    for i = 1: 1: n_sam
        new_r_mat(i, :) = new_r_mat(i, :) ./ sum_gmm_vec(i, 1);
    end
end

function pdf_vec = mgm (x_mat, mu_vec, sigma_mat)
    [n_dim, n_sam] = size(x_mat);
    
    x_mu_mat = zeros(n_dim, n_sam);
    for i = 1: 1: n_sam
        x_mu_mat(:, i) = x_mat(:, i) - mu_vec;
    end
    
    pdf_vec = exp((- 0.5) * sum(x_mu_mat .* (inv(sigma_mat) * x_mu_mat), 1)) / ...
              sqrt(((2 * pi) ^ n_dim) * det(sigma_mat)); %#ok<MINV>
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
