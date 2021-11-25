%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Get First N peaks of Angle Function
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

function ang_rad_x_vec = get_peak (abs_norm_vec, ang_rad_vec, n_src)
    [pks_vec, locs_x_vec] = findpeaks(abs_norm_vec);
    [~,       idx_x_vec]  = sort(pks_vec, 'descend');

    post_n_src = size(pks_vec, 1);
    if post_n_src > n_src
        post_n_src = n_src;
    end
    ang_rad_x_vec = ang_rad_vec(locs_x_vec(idx_x_vec(1: post_n_src, 1), 1), 1);

    n_lost_x = n_src - post_n_src;
    if n_lost_x > 0
        ang_rad_x_vec = [zeros(n_lost_x, 1); ang_rad_x_vec];
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%