% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;

% include
addpath('capture/');
addpath('miscellaneous/');
addpath('model/');
sam_path = 'output/record/d9-ch8-src3';

% time capture
is_pair     = true;
sim_max_tim = 10000e-3;

% real environment
n_src          = 3;
real_x_deg_vec = [ 48;  98; 152];
real_z_deg_vec = [ 86;  95;  95];

% common algorithm configuration
ang_rad_vec = (0: 1: 180).' .* pi / 180;

% test set
file_obj         = dir(fullfile(fullfile(pwd, sam_path), 'REC*.mat'));
sam_set_name_cel = {file_obj.name}.';
n_sam_set        = size(sam_set_name_cel, 1);

% overall configuration
n_mic     = 4;
n_mic_vec = [n_mic; 0; n_mic - 1];
c         = 343;
d_mic     = 0.09;
lamb      = 2 * d_mic;
cen_frq   = c / lamb;
err_val   = 1e-6;


% stft configuration
sam_frq      = 96e3;
sam_tim      = 1 / sam_frq;
win_fcn      = @blackman;
win_size     = fix(50e-3 / (1 / sam_frq)); % 100
n_fft        = (2 ^ 14);
hop_size     = ceil(win_size / (2 ^ 4));
stft_min_frq = 200;
stft_max_frq = cen_frq;

% number of points
n_pnt = fix(sim_max_tim * sam_frq);

% mic position
[~, pos_mat] = l_fwd_ster_fcn([n_mic; 0; n_mic - 1], d_mic, [ndef_ang, ndef_ang], lamb, false);

% test benchmarked
if is_pair == false
    real_x_deg_vec = sort(real_x_deg_vec, 'ascend');
    real_z_deg_vec = sort(real_z_deg_vec, 'ascend');
else
    [real_x_deg_vec, idx_vec] = sort(real_x_deg_vec, 'ascend');
    real_z_deg_vec            = real_z_deg_vec(idx_vec);
end
calc_x_rad_mat3 = zeros(n_src, 0, n_sam_set);
calc_z_rad_mat3 = zeros(n_src, 0, n_sam_set);
for i_sam_set = 1: 1: n_sam_set
    legend_name_cel = {};
    calc_x_rad_mat  = zeros(n_src, 0);
    calc_z_rad_mat  = zeros(n_src, 0);
    
    % load sample
    sam_set_name = sam_set_name_cel{i_sam_set, 1};
    load(fullfile(sam_path, sam_set_name));
    src_mat      = obj.wav_mat;
    src_mat      = src_mat(1: n_pnt, :);
    
    % sfft procedure
    [stft_mat3, frq_vec, tim_vec]   = stft_mat(src_mat, sam_frq, win_fcn, win_size, n_fft, hop_size);
    [post_stft_mat3, post_frq_vec]  = stft_hard_fitr(stft_mat3, frq_vec, stft_min_frq, stft_max_frq);
    post_stft_mat3                  = post_stft_mat3 ...
                                   ./ (win_size .* (sum(feval(win_fcn, win_size)) / win_size));
    
    if is_pair == false
        % doa procedure - no pair
        legend_name_cel{end + 1, 1} = 'IMUSIC';
        cd('doa/1983 - 12 - G. Su/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, ang_rad_vec);
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'TOFS';
        cd('doa/2007 - 10 - H. Yu/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, ang_rad_vec);
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'TOPS';
        cd('doa/2006 - 06 - Y. S. Yoon/');
        try
            [~, idx_ref] = min(abs(post_frq_vec - cen_frq));
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, idx_ref, ang_rad_vec);
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'Squared TOPS';
        cd('doa/2010 - 06 - K. Okane/');
        try
            [~, idx_ref] = min(abs(post_frq_vec - cen_frq));
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, idx_ref, ang_rad_vec);
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'WS-TOPS';
        cd('doa/2016 - 12 - H. Hayashi/');
        try
            n_post_frq = size(post_frq_vec, 1);
            idx_frq_ref_cvec = n_post_frq: 1: n_post_frq;
            tim_init = cputime;
            [calc_x_rad_vec, calc_z_rad_vec] = ...
                doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, idx_frq_ref_cvec, ang_rad_vec);
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'Proposed Method with MUSIC';
        cd('doa/now - B. Suksiri/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, "music-no-pair");
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'Proposed Method with ESPRIT';
        cd('doa/now - B. Suksiri/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = doa_frq_n_src(post_stft_mat3, post_frq_vec, n_mic, d_mic, lamb, n_src, cen_frq, "esprit-no-pair");
            calc_x_rad_vec = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec = sort(calc_z_rad_vec, 'ascend');
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');
    else
        % doa procedure - auto pair
        legend_name_cel{end + 1, 1} = '2D-IMUSIC';
        cd('doa/1983 - 12 - G. Su/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = abs_doa_frq_n_src(post_stft_mat3, post_frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec);
            [calc_x_rad_vec, idx_vec] = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec            = calc_z_rad_vec(idx_vec);
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = '2D-TOFS';
        cd('doa/2007 - 10 - H. Yu/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = abs_doa_frq_n_src(post_stft_mat3, post_frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec);
            [calc_x_rad_vec, idx_vec] = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec            = calc_z_rad_vec(idx_vec);
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');

        legend_name_cel{end + 1, 1} = 'Proposed Method with 2D-MUSIC';
        cd('doa/now - B. Suksiri/');
        try
            [calc_x_rad_vec, calc_z_rad_vec] = abs_doa_frq_n_src(post_stft_mat3, post_frq_vec, pos_mat, lamb, n_src, cen_frq, ang_rad_vec);
            [calc_x_rad_vec, idx_vec] = sort(calc_x_rad_vec, 'ascend');
            calc_z_rad_vec            = calc_z_rad_vec(idx_vec);
        catch
            calc_x_rad_vec = nan(n_src, 1);
            calc_z_rad_vec = nan(n_src, 1);
        end
        calc_x_rad_mat(:, end + 1) = calc_x_rad_vec;
        calc_z_rad_mat(:, end + 1) = calc_z_rad_vec;
        cd('../..');
    end
    
    n_test_set                                   = size(calc_x_rad_mat, 2);
    calc_x_rad_mat3(:, 1: n_test_set, i_sam_set) = calc_x_rad_mat;
    calc_z_rad_mat3(:, 1: n_test_set, i_sam_set) = calc_z_rad_mat;
end

to_str_fcn = @(in_str) num2str(in_str, '%.4f');

calc_x_deg_mat3 = calc_x_rad_mat3 .* 180 ./ pi;
calc_z_deg_mat3 = calc_z_rad_mat3 .* 180 ./ pi;

data_mean_x_deg_cel     = cell(n_test_set, n_src);
data_mean_z_deg_cel     = cell(n_test_set, n_src);
data_mean_deg_cel       = cell(n_test_set, 2 * n_src);
data_std_x_deg_cel      = cell(n_test_set, n_src);
data_std_z_deg_cel      = cell(n_test_set, n_src);
data_std_deg_cel        = cell(n_test_set, 2 * n_src);
data_rmse_x_deg_cel     = cell(n_test_set, n_src);
data_rmse_z_deg_cel     = cell(n_test_set, n_src);
data_rmse_deg_cel       = cell(n_test_set, 2 * n_src);
data_all_rmse_x_deg_cel = cell(n_test_set, 1);
data_all_rmse_z_deg_cel = cell(n_test_set, 1);
data_all_rmse_deg_cel   = cell(n_test_set, 1);
for i_test_set = 1: n_test_set
    all_err_x_deg_vec = [];
    all_err_z_deg_vec = [];
    
    for i_src = 1: n_src
        calc_x_deg_vec = squeeze(calc_x_deg_mat3(i_src, i_test_set, :));
        calc_z_deg_vec = squeeze(calc_z_deg_mat3(i_src, i_test_set, :));
        
        err_x_deg_vec = calc_x_deg_vec - real_x_deg_vec(i_src, 1);
        err_z_deg_vec = calc_z_deg_vec - real_z_deg_vec(i_src, 1);
        
        all_err_x_deg_vec = [all_err_x_deg_vec; err_x_deg_vec];
        all_err_z_deg_vec = [all_err_z_deg_vec; err_z_deg_vec];
        
        data_mean_x_deg_cel{i_test_set, i_src} = to_str_fcn(mean(calc_x_deg_vec));
        data_mean_z_deg_cel{i_test_set, i_src} = to_str_fcn(mean(calc_z_deg_vec));
        data_std_x_deg_cel {i_test_set, i_src} = to_str_fcn(std(calc_x_deg_vec));
        data_std_z_deg_cel {i_test_set, i_src} = to_str_fcn(std(calc_z_deg_vec));
        data_rmse_x_deg_cel{i_test_set, i_src} = to_str_fcn(mean(sqrt(err_x_deg_vec .^ 2)));
        data_rmse_z_deg_cel{i_test_set, i_src} = to_str_fcn(mean(sqrt(err_z_deg_vec .^ 2)));
        
        data_mean_deg_cel{i_test_set, (2 * i_src) - 1} = data_mean_x_deg_cel{i_test_set, i_src};
        data_mean_deg_cel{i_test_set, (2 * i_src) - 0} = data_mean_z_deg_cel{i_test_set, i_src};
        data_std_deg_cel {i_test_set, (2 * i_src) - 1} = data_std_x_deg_cel {i_test_set, i_src};
        data_std_deg_cel {i_test_set, (2 * i_src) - 0} = data_std_z_deg_cel {i_test_set, i_src};
        data_rmse_deg_cel{i_test_set, (2 * i_src) - 1} = data_rmse_x_deg_cel{i_test_set, i_src};
        data_rmse_deg_cel{i_test_set, (2 * i_src) - 0} = data_rmse_z_deg_cel{i_test_set, i_src};
    end
    data_all_rmse_x_deg_cel{i_test_set, 1} = mean(sqrt(all_err_x_deg_vec .^ 2));
    data_all_rmse_z_deg_cel{i_test_set, 1} = mean(sqrt(all_err_z_deg_vec .^ 2));
    data_all_rmse_deg_cel  {i_test_set, 1} = mean(sqrt([all_err_x_deg_vec; all_err_z_deg_vec] .^ 2));
end

writetable(cell2table([legend_name_cel, data_mean_deg_cel]),   strcat('doa_ap_3_real_mean',   '.xls'));
writetable(cell2table([legend_name_cel, data_mean_x_deg_cel]), strcat('doa_ap_3_real_mean_x', '.xls'));
writetable(cell2table([legend_name_cel, data_mean_z_deg_cel]), strcat('doa_ap_3_real_mean_z', '.xls'));
writetable(cell2table([legend_name_cel, data_std_deg_cel]),    strcat('doa_ap_3_real_std',    '.xls'));
writetable(cell2table([legend_name_cel, data_std_x_deg_cel]),  strcat('doa_ap_3_real_std_x',  '.xls'));
writetable(cell2table([legend_name_cel, data_std_z_deg_cel]),  strcat('doa_ap_3_real_std_z',  '.xls'));
writetable(cell2table([legend_name_cel, data_rmse_deg_cel,   data_all_rmse_deg_cel]),   strcat('doa_ap_3_real_rmse', '.xls'));
writetable(cell2table([legend_name_cel, data_rmse_x_deg_cel, data_all_rmse_x_deg_cel]), strcat('doa_ap_3_real_rmse_x', '.xls'));
writetable(cell2table([legend_name_cel, data_rmse_z_deg_cel, data_all_rmse_z_deg_cel]), strcat('doa_ap_3_real_rmse_z', '.xls'));
