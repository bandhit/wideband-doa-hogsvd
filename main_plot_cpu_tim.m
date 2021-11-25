% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;
version -blas
version -lapack

% include
addpath('miscellaneous/');

% test set
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

% test_set_name_cel{end + 1, 1} = 'out_now_b_suksiri_norm_mode_2.mat';
% legend_name_cel{end + 1, 1}   = 'Proposed Method with 2D-MUSIC';
% marker_str_cel{end + 1, 1}    = 'square';
% marker_size_cel{end + 1, 1}   = 10;

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

% path
out_path_str_cel = {};
n_mic_vec = [];
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch04';
n_mic_vec(end + 1, 1) = 4;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch05';
n_mic_vec(end + 1, 1) = 5;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch06';
n_mic_vec(end + 1, 1) = 6;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch07';
n_mic_vec(end + 1, 1) = 7;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch08';
n_mic_vec(end + 1, 1) = 8;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch09';
n_mic_vec(end + 1, 1) = 9;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch10';
n_mic_vec(end + 1, 1) = 10;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch11';
n_mic_vec(end + 1, 1) = 11;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch12';
n_mic_vec(end + 1, 1) = 12;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch13';
n_mic_vec(end + 1, 1) = 13;
out_path_str_cel{end + 1, 1} = 'output/computational complexity/d5-ch14';
n_mic_vec(end + 1, 1) = 14;

% plot
fig_obj   = figure;
file_name = 'cpu_tim_X';
x_val     = 21.0;
y_val     = 29.7;
subplot(3, 2, 1);
get_cpu_tim(out_path_str_cel, n_mic_vec, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel);
save_pdf(file_name, fig_obj, x_val, y_val);
pause(1);
close(fig_obj);

function get_cpu_tim(out_path_str_cel, n_mic_vec, test_set_name_cel, legend_name_cel, marker_str_cel, marker_size_cel)
    % plot configuration
    tim_min_val    = 0;
    tim_max_val    = 5;
    tim_step_val   = 1;
    n_sen_min_val  = 4;
    n_sen_max_val  = 14;
    n_sen_step_val = 2;
    
    % test benchmarked
    n_test_set = size(test_set_name_cel, 1);
    n_out_path = size(out_path_str_cel, 1);
    for i_test_set = 1: 1: n_test_set
        test_set_name = test_set_name_cel{i_test_set, 1};
        
        mean_diff_tim_vec = zeros(n_out_path, 1);
        std_diff_tim_vec  = zeros(n_out_path, 1);
        for i_out_path = 1: 1: n_out_path
            out_path_str = out_path_str_cel{i_out_path, 1};
            
            load(fullfile(out_path_str, test_set_name));

            mean_diff_tim_vec(i_out_path, 1) = mean(diff_tim_mat(:));
            std_diff_tim_vec(i_out_path, 1)  = std(diff_tim_mat(:));
        end
        
        legend_name = legend_name_cel{i_test_set, 1};
        marker_str  = marker_str_cel{i_test_set, 1};
        marker_size = marker_size_cel{i_test_set, 1};
        
        hold on;
        errorbar(n_mic_vec, mean_diff_tim_vec, std_diff_tim_vec, '-.', 'DisplayName', legend_name, 'Marker', marker_str, 'MarkerSize', marker_size, 'LineWidth', 0.75);
%         plot(n_mic_vec, mean_diff_tim_vec, '-.', 'DisplayName', legend_name, 'Marker', marker_str, 'MarkerSize', marker_size, 'LineWidth', 0.75);
        xlim([n_sen_min_val, n_sen_max_val]);
        xlabel('Number of Microphones, M', 'Color', 'k');
        ylabel('Runtime, Second', 'Color', 'k');
        grid on;
        set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'on');
        set(gca, 'XColor', 'k', 'XTick', n_sen_min_val: n_sen_step_val: n_sen_max_val, 'XTickLabel', string(n_sen_min_val: n_sen_step_val: n_sen_max_val));
        set(gca, 'FontName', 'Times New Roman');
        ylim([tim_min_val, tim_max_val]);
        set(gca, 'YColor', 'k', 'YTick', tim_min_val: tim_step_val: tim_max_val, 'YTickLabel', string(tim_min_val: tim_step_val: tim_max_val));
    end
%     legend('show', 'Location', 'southwest');
end
