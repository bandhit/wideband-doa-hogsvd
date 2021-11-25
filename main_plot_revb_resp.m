% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;

% include
addpath('miscellaneous/');
revb_path = 'reverberation';

fig_obj   = figure;
file_name = 'rir_resp';
x_val     = 24.0;
y_val     = 35.0;

wal_ref_mat = zeros(0, 6);

subplot(5, 2, 1);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.2-[15-15-5]-L12-S3-20190322T135417.mat', 1, 1);
subplot(5, 2, 2);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.3-[15-15-5]-L12-S3-20190320T223242.mat', 1, 1);
subplot(5, 2, 3);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.4-[15-15-5]-L12-S3-20190321T101148.mat', 1, 1);
subplot(5, 2, 4);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.5-[15-15-5]-L12-S3-20190321T155653.mat', 1, 1);
subplot(5, 2, 5);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.6-[15-15-5]-L12-S3-20190321T225903.mat', 1, 1);
subplot(5, 2, 6);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.7-[15-15-5]-L12-S3-20190322T054003.mat', 1, 1);
subplot(5, 2, 7);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.8-[15-15-5]-L12-S3-20190323T174940.mat', 1, 1);
subplot(5, 2, 8);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-0.9-[15-15-5]-L12-S3-20190323T161401.mat', 1, 1);
subplot(5, 2, 9);
wal_ref_mat(end + 1, 1: 6) = plot_resp(revb_path, 'T60-1-[15-15-5]-L12-S3-20190325T120146.mat', 1, 1);

% save_pdf(file_name, fig_obj, x_val, y_val);
% pause(1);
% close(fig_obj);

writetable(array2table(wal_ref_mat), strcat('wal_ref_rir', '.xls'))

function wal_ref_cvec = plot_resp(revb_path, revb_file_name_string, idx_src, idx_sum_mic)
    load(fullfile(revb_path, revb_file_name_string), 'rir_cel', 'sam_frq', 'rt_val', 'wal_ref_cvec', 'room_dim_cvec', 'alpha_cvec', 'c');
    rir_vec      = rir_cel{idx_src, idx_sum_mic};
%     log_rir_vec  = 10 .* log10(abs(rir_vec));
%     log_rir_vec  = log_rir_vec - (max(log_rir_vec) / 2);
    log_rir_vec  = 10 .* log10(abs(rir_vec)) + 5;
    n_sam        = size(rir_vec, 1);
    tim_vec      = (0: 1: n_sam - 1).' ./ sam_frq;
    offset_rt    = -60;
    offset_rir   = -80;
    offset_tim   = 0.25 * rt_val;
    
    plot(tim_vec, log_rir_vec, 'Color', 'b');
    xlim([0, rt_val + offset_tim]);
    ylim([offset_rir, 0]);
    grid on;
    grid minor;
    hold on;
    plot([rt_val; rt_val], [0; offset_rir], 'Color', 'k');
    plot([0; rt_val + offset_tim], [offset_rt; offset_rt], 'Color', 'k');
    plot(tim_vec, (offset_rt / sqrt(rt_val)) * sqrt(tim_vec), 'Color', 'r');
    text(rt_val, 0.33 * offset_rir, strcat('\leftarrow ', string(rt_val), ' s'));
    text(1.06 * rt_val, 0.6 * offset_rir, 'RT60');
    xlabel(['Time (Second)', newline, '(?)']);
    ylabel('RIR (Decibel)');
    
%     % alpha_cvec is absorption coefficients, so wal_ref_cvec is reflection coefficient
%     load(fullfile(revb_path, revb_file_name_string), 'wal_ref_cvec', 'room_dim_cvec', 'alpha_cvec', 'c');
%     sab_rt_val = sab_rt_fcn(c, room_dim_cvec.', alpha_cvec.');
%     eyr_rt_val = eyr_rt_fcn(c, room_dim_cvec.', alpha_cvec.');
end
