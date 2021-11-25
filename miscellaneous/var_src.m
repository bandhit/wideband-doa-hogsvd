%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Variance of Sources Calculation
%
% Description : Variance of Sources Calculation
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 9 September 2017, Bandhit Suksiri,
%               Updated: 9 September 2017, Bandhit Suksiri.
%
% Copyright 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function meas_s_mat = var_src (stft_mat3)
    [n_frq, n_tim, n_src] = size(stft_mat3);
    meas_s_mat            = zeros(n_src, n_frq);
    for i_frq = 1: 1: n_frq
        stft_mat = squeeze(stft_mat3(i_frq, :, :)).';
        
        corr_s_mat = zeros(n_src, n_src);
        for i_tim = 1: 1: n_tim
            corr_s_mat = corr_s_mat + (stft_mat(:, i_tim) * stft_mat(:, i_tim)');
        end
        corr_s_mat = corr_s_mat / n_tim;
        meas_s_mat(:, i_frq) = real(diag(corr_s_mat));
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%