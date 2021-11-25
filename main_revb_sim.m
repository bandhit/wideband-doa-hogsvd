% reset workspace
path(pathdef);
close all;
clear playsnd;
clear;
clc;

% include
addpath('simulation/');
addpath('model/');
addpath('miscellaneous/');
addpath('third party/Fast simulation of acoustic room impulse responses/');

fprintf('>> RIR 0.2 RT60\n');
revb_run(0.2, [13, 13, 5], 120, 60, 'T60', 0.065);
fprintf('>> RIR 0.3 RT60\n');
revb_run(0.3, [15, 15, 5], 90, 30, 'T60', 0.035);
fprintf('>> RIR 0.4 RT60\n');
revb_run(0.4, [15, 15, 5], 90, 30, 'T60', 0.035);
fprintf('>> RIR 0.5 RT60\n');
revb_run(0.5, [15, 15, 5], 90, 30, 'T60', 0.035);
fprintf('>> RIR 0.6 RT60\n');
revb_run(0.6, [15, 15, 5], 90, 30, 'T60', 0.035);
fprintf('>> RIR 0.7 RT60\n');
revb_run(0.7, [15, 15, 5], 90, 30, 'T60', 0.035);
fprintf('>> RIR 0.8 RT60\n');
revb_run(0.8, [15, 15, 5], 90, 30, 'T60', 0.035);

function revb_run (rt_val, room_dim_cvec, delta_db_val, diffu_db_val, rt_type_str, con_std_val)
    % include
    revb_path = 'reverberation';

    % overall configuration
    n_mic        = 12;
    n_mic_vec    = [n_mic; 0; n_mic - 1];
    c            = 343;
    d_mic        = 0.05; % 0.09
    lamb         = 2 * d_mic;
    cen_frq      = c / lamb;
    err_val      = 1e-12;
    dis_terb_val = 1;

    % source position configuration
    az_deg_vec  = [60; 45; 30];
    el_deg_vec  = [30; 45; 60];
    amp_src_vec = [ 3;  3;  3];

    % room reverberation configuration
    sam_frq     = 192e3;
    is_orth_pnt = true;

    % simulation
    sim_obj = wav_sim_3d_def();
    for info_vec = [az_deg_vec, el_deg_vec, amp_src_vec].'
        az_deg  = info_vec(1, 1);
        el_deg  = info_vec(2, 1);
        amp_src = info_vec(3, 1);
        sim_obj.add_src_ang(az_deg - err_val, el_deg - err_val, amp_src);
    end
    [~, pos_mat] = l_fwd_ster_fcn(n_mic_vec, d_mic, [ndef_ang, ndef_ang], lamb, false);
    pos_mat = pos_mat + dis_terb_val;
    for pos_vec = pos_mat.'
        sim_obj.add_pnt(pos_vec);
    end
    sim_obj.prep();

    % wall reflection calculation
    is_pass = 0;
    while is_pass == 0
        wight_cvec = (rand(1, 6) .* 0.9 + 0.1) .^ rand;
        [alpha_cvec, is_pass] = ISM_AbsCoeff(rt_type_str, rt_val, room_dim_cvec, wight_cvec, 'lehmannjohansson'); % sabine
        if is_pass == 1
            surf_x       = room_dim_cvec(1, 2) * room_dim_cvec(1, 3);
            surf_y       = room_dim_cvec(1, 1) * room_dim_cvec(1, 3);
            surf_z       = room_dim_cvec(1, 1) * room_dim_cvec(1, 2);
            avec         = [surf_x * alpha_cvec(1, 1), ...
                            surf_x * alpha_cvec(1, 2), ...
                            surf_y * alpha_cvec(1, 3), ...
                            surf_y * alpha_cvec(1, 4), ...
                            surf_z * alpha_cvec(1, 5), ...
                            surf_z * alpha_cvec(1, 6)];
            std_val      = std(avec / (2 * (surf_x + surf_y + surf_z)));
            if (std_val > con_std_val) % 0.01 for 200 ms, 0.035 for >200 ms
                is_pass = 0;
            end
        end
    end
    wal_ref_cvec = sqrt(1 - alpha_cvec); % alpha_cvec is absorption coefficients, so wal_ref_cvec is reflection coefficients
    
    % old-style RT60 measurement
    sab_rt_val = sab_rt_fcn(c, room_dim_cvec.', alpha_cvec.');
    eyr_rt_val = eyr_rt_fcn(c, room_dim_cvec.', alpha_cvec.');
    
    % reverberation
    orth_pnt_mat3 = sim_obj.orth_pnt_mat3;
    cen_pos_mat   = sim_obj.cen_pos_mat;
    [~, n_src, sum_mic] = size(orth_pnt_mat3);
    rir_cel = cell(n_src, sum_mic);
    for i_sum_mic = 1: 1: sum_mic
        mic_pos_cvec = pos_mat(i_sum_mic, :);
        parfor i_src = 1: 1: n_src
            if is_orth_pnt == true
                src_pos_cvec = orth_pnt_mat3(:, i_src, i_sum_mic).';
            else
                src_pos_cvec = cen_pos_mat(:, i_src).';
            end
            src_pos_cvec = src_pos_cvec + err_val;
            rir_vec = fast_ISM_RoomResp(sam_frq, wal_ref_cvec, rt_type_str, rt_val, src_pos_cvec, mic_pos_cvec, room_dim_cvec, 'c', c, 'Delta_dB', delta_db_val, 'Diffuse_dB', diffu_db_val);
            rir_cel{i_src, i_sum_mic} = rir_vec;
            
%             % plot test
%             rir_vec = abs(rir_vec);
%             rir_vec = rir_vec ./ max(rir_vec);
%             figure;
%             plot((1: 1: size(rir_vec, 1)) / sam_frq, 20 * log10(rir_vec));
%             figure;
%             plot((1: 1: size(rir_vec, 1)) / sam_frq, rir_vec);
            
            fprintf('>> %6.2f\n', (((i_sum_mic - 1) * n_src) + i_src) / (sum_mic * n_src) * 100);
        end
    end

    clear sim_obj;
    file_name_string = char(strcat(rt_type_str, '-', ...
                              string(rt_val), '-', ...
                              '[', string(room_dim_cvec(1, 1)), '-', string(room_dim_cvec(1, 2)), '-', string(room_dim_cvec(1, 3)), ']', '-', ...
                              'L', string(n_mic), '-', ...
                              'S', string(n_src), '-', ...
                              datestr(now,'yyyymmddTHHMMSS'), ...
                              '.mat'));
    save(fullfile(revb_path, file_name_string));
end
