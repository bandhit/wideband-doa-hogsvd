% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;

% include
addpath('miscellaneous/');

% test set 1
if true
    test_set_name_cel             = {};
    legend_name_cel               = {};
    marker_str_cel                = {};
    marker_size_cel               = {};

    test_set_name_cel{end + 1, 1} = 'out_now_b_suksiri_norm_mode_0.mat';
    legend_name_cel{end + 1, 1}   = 'Proposed Method with MUSIC';
    marker_str_cel{end + 1, 1}    = 'o';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_now_b_suksiri_norm_mode_1.mat';
    legend_name_cel{end + 1, 1}   = 'Proposed Method with ESPRIT';
    marker_str_cel{end + 1, 1}    = 'x';
    marker_size_cel{end + 1, 1}   = 10;
    
    test_set_name_cel{end + 1, 1} = 'out_1983_12_g_su.mat';
    legend_name_cel{end + 1, 1}   = 'IMUSIC';
    marker_str_cel{end + 1, 1}    = 'diamond';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_2007_10_h_yu.mat';
    legend_name_cel{end + 1, 1}   = 'TOFS';
    marker_str_cel{end + 1, 1}    = '+';
    marker_size_cel{end + 1, 1}   = 10;

    test_set_name_cel{end + 1, 1} = 'out_2006_06_y_s_yoon.mat';
    legend_name_cel{end + 1, 1}   = 'TOPS';
    marker_str_cel{end + 1, 1}    = '^';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_2010_06_k_okane.mat';
    legend_name_cel{end + 1, 1}   = 'Squared TOPS';
    marker_str_cel{end + 1, 1}    = 'hexagram';
    marker_size_cel{end + 1, 1}   = 9;

    test_set_name_cel{end + 1, 1} = 'out_2016_12_h_hayashi.mat';
    legend_name_cel{end + 1, 1}   = 'WS-TOPS';
    marker_str_cel{end + 1, 1}    = '>';
    marker_size_cel{end + 1, 1}   = 7;

    % plot
    fig_obj = figure;
    file_name = 'rmse_n_mic_no_pair';
    x_val     = 29.7;
    y_val     = 21.0 - 4;
    subplot(2, 3, 1);
    plot_rmse('output/number of mic and distance/no pair/d9-ch04', true,  30, -1, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(2, 3, 2);
    plot_rmse('output/number of mic and distance/no pair/d9-ch08', true,  30, -1, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(2, 3, 3);
    plot_rmse('output/number of mic and distance/no pair/d9-ch12', true,  30, -1, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    
    subplot(2, 3, 4);
    plot_rmse('output/number of mic and distance/no pair/d9-ch04', false, 30, -2, true,  test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(2, 3, 5);
    plot_rmse('output/number of mic and distance/no pair/d9-ch08', false, 30, -2, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(2, 3, 6);
    plot_rmse('output/number of mic and distance/no pair/d9-ch12', false, 30, -2, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    save_pdf(file_name, fig_obj, x_val, y_val);
    pause(1);
    close(fig_obj);
end

% test set 2
if true
    test_set_name_cel             = {};
    legend_name_cel               = {};
    marker_str_cel                = {};
    marker_size_cel               = {};

    test_set_name_cel{end + 1, 1} = 'out_now_b_suksiri_norm_mode_2.mat';
    legend_name_cel{end + 1, 1}   = 'Proposed Method with 2D-MUSIC';
    marker_str_cel{end + 1, 1}    = 'square';
    marker_size_cel{end + 1, 1}   = 10;

    test_set_name_cel{end + 1, 1} = 'out_1983_12_g_su.mat';
    legend_name_cel{end + 1, 1}   = '2D-IMUSIC';
    marker_str_cel{end + 1, 1}    = 'diamond';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_2007_10_h_yu.mat';
    legend_name_cel{end + 1, 1}   = '2D-TOFS';
    marker_str_cel{end + 1, 1}    = '+';
    marker_size_cel{end + 1, 1}   = 10;

    % plot
    fig_obj = figure;
    file_name = 'rmse_n_mic_auto_pair';
    x_val     = 21.0;
    y_val     = 29.7;
    subplot(3, 2, 1);
    plot_rmse('output/number of mic and distance/auto pair/d5-ch04', true,  40, -1, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(3, 2, 3);
    plot_rmse('output/number of mic and distance/auto pair/d5-ch04', false, 40, -2, true,  test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    save_pdf(file_name, fig_obj, x_val, y_val);
    pause(1);
    close(fig_obj);
end

% test set 3
if true
    test_set_name_cel             = {};
    legend_name_cel               = {};
    marker_str_cel                = {};
    marker_size_cel               = {};

    test_set_name_cel{end + 1, 1} = 'out_now_b_suksiri_norm_mode_0.mat';
    legend_name_cel{end + 1, 1}   = 'Proposed Method with MUSIC';
    marker_str_cel{end + 1, 1}    = 'o';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_now_b_suksiri_norm_mode_1.mat';
    legend_name_cel{end + 1, 1}   = 'Proposed Method with ESPRIT';
    marker_str_cel{end + 1, 1}    = 'x';
    marker_size_cel{end + 1, 1}   = 10;
    
    test_set_name_cel{end + 1, 1} = 'out_1983_12_g_su.mat';
    legend_name_cel{end + 1, 1}   = 'IMUSIC';
    marker_str_cel{end + 1, 1}    = 'diamond';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_2007_10_h_yu.mat';
    legend_name_cel{end + 1, 1}   = 'TOFS';
    marker_str_cel{end + 1, 1}    = '+';
    marker_size_cel{end + 1, 1}   = 10;

    test_set_name_cel{end + 1, 1} = 'out_2006_06_y_s_yoon.mat';
    legend_name_cel{end + 1, 1}   = 'TOPS';
    marker_str_cel{end + 1, 1}    = '^';
    marker_size_cel{end + 1, 1}   = 8;

    test_set_name_cel{end + 1, 1} = 'out_2010_06_k_okane.mat';
    legend_name_cel{end + 1, 1}   = 'Squared TOPS';
    marker_str_cel{end + 1, 1}    = 'hexagram';
    marker_size_cel{end + 1, 1}   = 9;

    test_set_name_cel{end + 1, 1} = 'out_2016_12_h_hayashi.mat';
    legend_name_cel{end + 1, 1}   = 'WS-TOPS';
    marker_str_cel{end + 1, 1}    = '>';
    marker_size_cel{end + 1, 1}   = 7;

    % plot
    fig_obj = figure;
    file_name = 'rmse_src';
    x_val     = 21.0;
    y_val     = 29.7;
    subplot(3, 2, 1);
    plot_rmse('output/verious source type/sine/d5-ch06',  true,  40, -1, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(3, 2, 2);
    plot_rmse('output/verious source type/human/d5-ch06', true,  40, -1, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(3, 2, 3);
    plot_rmse('output/verious source type/sine/d5-ch06',  false, 40, -2, true,  test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    subplot(3, 2, 4);
    plot_rmse('output/verious source type/human/d5-ch06', false, 40, -2, false, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
    save_pdf(file_name, fig_obj, x_val, y_val);
    pause(1);
    close(fig_obj);
end

function plot_rmse(out_path, is_rmse, snr_max_val, min_val, is_lgnd, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel)
    % plot configuration
    snr_min_val   = -10;
%     snr_max_val   = 40;
    snr_step_val  = 10;
    std_min_val   = min_val;
%     std_min_val   = -3;
    std_max_val   = 2;
    std_step_val  = 1;
    rmse_min_val  = min_val;
%     rmse_min_val  = -1;
    rmse_max_val  = 2;
    rmse_step_val = 1;
    skip_ratio    = 1;
    
    % test benchmarked
    n_test_set = size(test_set_name_cel, 1);
    for i_test_set = 1: 1: n_test_set
        test_set_name = test_set_name_cel{i_test_set, 1};

        load(fullfile(out_path, test_set_name));

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

        legend_name = legend_name_cel{i_test_set, 1};
        marker_str  = marker_str_cel{i_test_set, 1};
        marker_size = 0.75 * marker_size_cel{i_test_set, 1};

        hold on;
        if is_rmse == true
            plot(snr_vec([1: skip_ratio: end, end], 1), rmse_sam_vec([1: skip_ratio: end, end], 1), '-.', 'DisplayName', legend_name, 'Marker', marker_str, 'MarkerSize', marker_size, 'LineWidth', 0.75);
            set(gca, 'YScale', 'log');
            xlim([snr_min_val, snr_max_val]);
            ylim(10 .^ [rmse_min_val, rmse_max_val]);
            ylabel('RMSE, Degree', 'Color', 'k');
            xlabel({'SNR, dB'; '(?)'}, 'Color', 'k');
            grid on;
            set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'on');
            set(gca, 'XColor', 'k', 'XTick', snr_min_val: snr_step_val: snr_max_val, 'XTickLabel', string(snr_min_val: snr_step_val: snr_max_val));
            set(gca, 'YColor', 'k', 'YTick', 10 .^ (rmse_min_val: rmse_step_val: rmse_max_val), 'YTickLabel', strcat(strcat(string('10^{'), string(rmse_min_val: rmse_step_val: rmse_max_val)), string('}')));
            set(gca, 'FontName', 'Times New Roman');
        else
            plot(snr_vec([1: skip_ratio: end, end], 1), std_sam_vec([1: skip_ratio: end, end], 1), '-.', 'DisplayName', legend_name, 'Marker', marker_str, 'MarkerSize', marker_size, 'LineWidth', 0.75);
            set(gca, 'YScale', 'log');
            xlim([snr_min_val, snr_max_val]);
            ylim(10 .^ [std_min_val, std_max_val]);
            ylabel('SD, Degree', 'Color', 'k');
            xlabel({'SNR, dB'; '(?)'}, 'Color', 'k');
            grid on;
            set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'on');
            set(gca, 'XColor', 'k', 'XTick', snr_min_val: snr_step_val: snr_max_val, 'XTickLabel', string(snr_min_val: snr_step_val: snr_max_val));
            set(gca, 'YColor', 'k', 'YTick', 10 .^ (std_min_val: std_step_val: std_max_val), 'YTickLabel', strcat(strcat(string('10^{'), string(std_min_val: std_step_val: std_max_val)), string('}')));
            set(gca, 'FontName', 'Times New Roman');
        end
    end
    if is_lgnd == true
        legend('show', 'Location', 'southwest');
    end
end
