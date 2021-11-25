% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;

% include
addpath('miscellaneous/');

% reverberation
rt_vec   = [0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0];
out_path = 'output/reverberation/d5-ch12-rt60';

% plot
fig_obj   = figure;
file_name = 'rmse_revb_t60';
x_val     = 21.0 + 4;
y_val     = 29.7;

subplot(7, 2, 1);
plot_revb(out_path, true,  'out_1983_12_g_su.mat', rt_vec, 'IMUSIC');
subplot(7, 2, 2);
plot_revb(out_path, false, 'out_1983_12_g_su.mat', rt_vec, 'IMUSIC');    

subplot(7, 2, 3);
plot_revb(out_path, true,  'out_2007_10_h_yu.mat', rt_vec, 'TOFS');
subplot(7, 2, 4);
plot_revb(out_path, false, 'out_2007_10_h_yu.mat', rt_vec, 'TOFS');

subplot(7, 2, 5);
plot_revb(out_path, true,  'out_2006_06_y_s_yoon.mat', rt_vec, 'TOPS');
subplot(7, 2, 6);
plot_revb(out_path, false, 'out_2006_06_y_s_yoon.mat', rt_vec, 'TOPS');

subplot(7, 2, 7);
plot_revb(out_path, true,  'out_2010_06_k_okane.mat', rt_vec, 'Squared TOPS');
subplot(7, 2, 8);
plot_revb(out_path, false, 'out_2010_06_k_okane.mat', rt_vec, 'Squared TOPS');

subplot(7, 2, 9);
plot_revb(out_path, true,  'out_2016_12_h_hayashi.mat', rt_vec, 'WS-TOPS');
subplot(7, 2, 10);
plot_revb(out_path, false, 'out_2016_12_h_hayashi.mat', rt_vec, 'WS-TOPS');

subplot(7, 2, 11);
plot_revb(out_path, true,  'out_now_b_suksiri_norm_mode_0.mat', rt_vec, 'Proposed Method with MUSIC');
subplot(7, 2, 12);
plot_revb(out_path, false, 'out_now_b_suksiri_norm_mode_0.mat', rt_vec, 'Proposed Method with MUSIC');

subplot(7, 2, 13);
plot_revb(out_path, true,  'out_now_b_suksiri_norm_mode_1.mat', rt_vec, 'Proposed Method with ESPRIT');
xlabel({'SNR, dB'; '(?)'}, 'Color', 'k');
subplot(7, 2, 14);
plot_revb(out_path, false, 'out_now_b_suksiri_norm_mode_1.mat', rt_vec, 'Proposed Method with ESPRIT');
xlabel({'SNR, dB'; '(?)'}, 'Color', 'k');

% save_pdf(file_name, fig_obj, x_val, y_val);
% pause(1);
% close(fig_obj);
% 
% fig_obj_axis   = figure;
% file_name_axis = 'rmse_revb_rt60_axis';
% subplot(7, 2, 13);
% plot_revb(out_path, true, 'out_now_b_suksiri_norm_mode_0.mat', rt_vec, '');
% rmse_min_val  = 0;
% rmse_max_val  = 2;
% rmse_step_val = 1;
% bar = colorbar('Location', 'northoutside', 'YTick', rmse_min_val: rmse_step_val: rmse_max_val, 'YTickLabel', cellstr(strcat(strcat(string('10^{'), string(rmse_min_val: rmse_step_val: rmse_max_val)), string('}'))));
% ylabel(bar, 'Magnitude');
% subplot(7, 2, 14);
% plot_revb(out_path, false, 'out_now_b_suksiri_norm_mode_0.mat', rt_vec, '');
% std_min_val   = -1;
% std_max_val   = 2;
% std_step_val  = 1;
% bar = colorbar('Location', 'northoutside', 'YTick', std_min_val: std_step_val: std_max_val, 'YTickLabel', cellstr(strcat(strcat(string('10^{'), string(std_min_val: std_step_val: std_max_val)), string('}'))));
% ylabel(bar, 'Magnitude');
% save_pdf(file_name_axis, fig_obj_axis, x_val, y_val);
% pause(1);
% close(fig_obj_axis);

function plot_revb(out_path, is_rmse, test_set_name, rt_vec, legend_name)
    % plot configuration
    rt_min_val   = min(rt_vec);
    rt_max_val   = max(rt_vec);
    rt_step_val  = 0.2;
    snr_min_val  = -10;
    snr_max_val  = 40;
    snr_step_val = 10;
    rmse_min_val = 0;
    rmse_max_val = 2;
    std_min_val  = -1;
    std_max_val  = 2;
    
    % test benchmarked
    tmp_rt_str    = num2str(rt_vec(1, 1), '%.1f');
    load(fullfile(out_path, tmp_rt_str, test_set_name));
    [~, n_snr, ~] = size(calc_x_deg_mat3);
    n_rt          = size(rt_vec, 1);
    rmse_sam_mat  = zeros(n_snr, n_rt);
    std_sam_mat   = zeros(n_snr, n_rt);
    out_snr_vec   = snr_vec;
    for i_rt = 1: 1: n_rt
        tmp_rt_str = num2str(rt_vec(i_rt, 1), '%.1f');
        load(fullfile(out_path, tmp_rt_str, test_set_name));

        [n_src, n_snr, n_sam] = size(calc_x_deg_mat3);
        meas_x_deg_vec        = sort(meas_x_deg_vec, 'ascend');
        meas_z_deg_vec        = sort(meas_z_deg_vec, 'ascend');
        for i_snr = 1: 1: n_snr
            for i_sam = 1: 1: n_sam
                calc_x_deg_mat3(:, i_snr, i_sam) = real(sort(calc_x_deg_mat3(:, i_snr, i_sam), 'ascend'));
                calc_z_deg_mat3(:, i_snr, i_sam) = real(sort(calc_z_deg_mat3(:, i_snr, i_sam), 'ascend'));
            end
        end

        std_sam_vec  = zeros(n_snr, 1);
        rmse_sam_vec = zeros(n_snr, 1);
        for i_snr = 1: 1: n_snr
            dev_vec = [];
            err_vec = [];
            for i_src = 1: 1: n_src
                meas_x_deg     = meas_x_deg_vec(i_src, 1);
                meas_z_deg     = meas_z_deg_vec(i_src, 1);
                calc_x_deg_vec = squeeze(real(calc_x_deg_mat3(i_src, i_snr, :)));
                calc_z_deg_vec = squeeze(real(calc_z_deg_mat3(i_src, i_snr, :)));
                mean_sam_x     = mean(calc_x_deg_vec);
                mean_sam_z     = mean(calc_z_deg_vec);
                dev_vec        = [dev_vec; calc_x_deg_vec - mean_sam_x; calc_z_deg_vec - mean_sam_z];
                err_vec        = [err_vec; calc_x_deg_vec - meas_x_deg; calc_z_deg_vec - meas_z_deg];
            end
            std_sam  = mean(sqrt(dev_vec .^ 2));
            rmse_sam = mean(sqrt(err_vec .^ 2));
            if std_sam < eps
                std_sam = eps;
            end
            if rmse_sam < eps
                rmse_sam = eps;
            end
            std_sam_vec(i_snr, 1)  = std_sam;
            rmse_sam_vec(i_snr, 1) = rmse_sam;
        end
        
        std_sam_mat(:, i_rt)  = std_sam_vec;
        rmse_sam_mat(:, i_rt) = rmse_sam_vec;
    end
    
    if is_rmse == true
        obj = pcolor(out_snr_vec, rt_vec, log10(rmse_sam_mat.'));
        xlim([snr_min_val, snr_max_val]);
        ylim([rt_min_val,  rt_max_val]);
        set(gca, 'XColor', 'k', 'XTick', snr_min_val: snr_step_val: snr_max_val, 'XTickLabel', string(snr_min_val: snr_step_val: snr_max_val));
        set(gca, 'YColor', 'k', 'YTick', rt_min_val:  rt_step_val:  rt_max_val,  'YTickLabel', string(rt_min_val:  rt_step_val:  rt_max_val));
        caxis([rmse_min_val, rmse_max_val]);
    else
        obj = pcolor(out_snr_vec, rt_vec, log10(std_sam_mat.'));
        xlim([snr_min_val, snr_max_val]);
        ylim([rt_min_val,  rt_max_val]);
        set(gca, 'XColor', 'k', 'XTick', snr_min_val: snr_step_val: snr_max_val, 'XTickLabel', string(snr_min_val: snr_step_val: snr_max_val));
        set(gca, 'YColor', 'k', 'YTick', rt_min_val:  rt_step_val:  rt_max_val,  'YTickLabel', string(rt_min_val:  rt_step_val:  rt_max_val));
        caxis([std_min_val, std_max_val]);
    end
    set(gca, 'layer', 'top');
    grid on;
    shading(gca, 'faceted');
    colormap(jet(512));
%     set(obj, 'EdgeColor', 'none');
    set(obj, 'EdgeColor', [1, 1, 1]);
    set(gca, 'FontName', 'Times New Roman');
    
    ylabel({legend_name; 'RT60, Second'}, 'Color', 'k');
end
