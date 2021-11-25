% overall configuration
n_mic     = 4;
n_mic_vec = [n_mic; 0; n_mic - 1];
c         = 343;
d_mic     = 0.05; % 0.09
lamb      = 2 * d_mic;
cen_frq   = c / lamb;
err_val   = 1e-6;

% stft configuration
sam_frq      = 192e3; % 48k
sam_tim      = 1 / sam_frq;
win_fcn      = @blackman;
win_size     = fix(50e-3 / (1 / sam_frq)); % 100
n_fft        = (2 ^ 14);
hop_size     = ceil(win_size / (2 ^ 4));
% stft_min_frq = cen_frq - 200;
% stft_max_frq = cen_frq;
% stft_min_frq = 200;
% stft_max_frq = 3400;
% stft_min_frq = 200;
% stft_max_frq = 4000;
stft_min_frq = 200;
stft_max_frq = cen_frq;
% stft_min_frq = 700;
% stft_max_frq = 800;

% common algorithm configuration
ang_rad_vec = (0: 1: 180).' .* pi / 180;
is_pair = false;

% source position configuration
az_deg_vec  = [30; 45; 60];
el_deg_vec  = [60; 45; 30];
amp_src_vec = [ 1;  1;  1];
% az_deg_vec  = [80; 20];
% el_deg_vec  = [20; 80];
% amp_src_vec = [ 1;  1];

% source wave configuration
addpath('sample/fukulab-001/192k/'); % 48k
addpath('sample/4409__pinkyfinger__piano-notes-1-octave/frequency conversion/48k/');
sim_max_tim = 1000e-3;
ex_1_vec    = audioread('001-001-002-001-002.wav');
ex_2_vec    = audioread('003-001-002-006-001.wav');
% ex_3_vec    = audioread('004-001-002-002-001.wav');
ex_3_vec    = audioread('002-001-002-005-001.wav');
ex_1_vec    = ex_1_vec / min(abs([min(ex_1_vec), max(ex_1_vec)]));
ex_2_vec    = ex_2_vec / min(abs([min(ex_2_vec), max(ex_2_vec)]));
ex_3_vec    = ex_3_vec / min(abs([min(ex_3_vec), max(ex_3_vec)]));
raw_tim_vec = (0: sam_tim: sim_max_tim / 2)';
ex_4_vec    = rand(size(raw_tim_vec, 1), 1) .* sin(2 .* pi .* 0.78e3 .* raw_tim_vec);
ex_5_vec    = rand(size(raw_tim_vec, 1), 1) .* sin(2 .* pi .* 0.53e3 .* raw_tim_vec);
ex_6_vec    = rand(size(raw_tim_vec, 1), 1) .* sin(2 .* pi .* 0.44e3 .* raw_tim_vec);
% ex_4_vec    = sin(2 .* pi .* 3.20e3 .* raw_tim_vec);
% ex_5_vec    = sin(2 .* pi .* 1.70e3 .* raw_tim_vec);
% ex_6_vec    = sin(2 .* pi .* 0.85e3 .* raw_tim_vec);
ex_7_vec    = audioread('68448__pinkyfinger__piano-g.wav'); % 780 Hz (G5)
ex_8_vec    = audioread('68441__pinkyfinger__piano-c.wav'); % 530 Hz (C5, Tenor C)
ex_9_vec    = audioread('68437__pinkyfinger__piano-a.wav'); % 440 Hz (A4, A440)
ex_7_vec    = ex_7_vec / min(abs([min(ex_7_vec), max(ex_7_vec)]));
ex_8_vec    = ex_8_vec / min(abs([min(ex_8_vec), max(ex_8_vec)]));
ex_9_vec    = ex_9_vec / min(abs([min(ex_9_vec), max(ex_9_vec)]));
wav_cel     = {ex_1_vec, ex_2_vec, ex_3_vec};
% wav_cel     = {ex_4_vec, ex_5_vec, ex_6_vec};
% wav_cel     = {ex_7_vec, ex_8_vec, ex_9_vec};
% wav_cel     = {ex_1_vec, ex_2_vec};
% wav_cel     = {ex_2_vec, ex_3_vec};
% wav_cel     = {ex_7_vec, ex_8_vec};

% room reverberation configuration
revb_file_name_string = NaN;
% revb_file_name_string = 'T60-0.2-[15-15-5]-L12-S3-20190322T135417.mat';
% revb_file_name_string = 'T60-0.3-[15-15-5]-L12-S3-20190320T223242.mat';
% revb_file_name_string = 'T60-0.4-[15-15-5]-L12-S3-20190321T101148.mat';
% revb_file_name_string = 'T60-0.5-[15-15-5]-L12-S3-20190321T155653.mat';
% revb_file_name_string = 'T60-0.6-[15-15-5]-L12-S3-20190321T225903.mat';
% revb_file_name_string = 'T60-0.7-[15-15-5]-L12-S3-20190322T054003.mat';
% revb_file_name_string = 'T60-0.8-[15-15-5]-L12-S3-20190323T174940.mat';
% revb_file_name_string = 'T60-0.9-[15-15-5]-L12-S3-20190323T161401.mat';
% revb_file_name_string = 'T60-1-[15-15-5]-L12-S3-20190325T120146.mat';

% test procedure configuration
% n_sam   = 100;
% snr_vec = (-10: 5: 40).';
n_sam   = 1;
snr_vec = [10; 20; 30; 40; 50];
% n_sam   = 1;
% snr_vec = [40];
% n_sam   = 100;
% snr_vec = [100];

% warning ('off', 'all');
