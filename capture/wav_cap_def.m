%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Class of Audio Interface Capture
%
% Description : Audio Interface Capture (OS: Windows)
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 16 May 2018, Bandhit Suksiri,
%               Updated: 17 May 2018, Bandhit Suksiri.
%
% Copyright 2016 - 2018,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef wav_cap_def < handle
    properties (Constant = true)
        DIVR_CONF_DIC_SND_STR = 'DirectSound';
        DIVR_CONF_ASIO_STR    = 'ASIO';
        DIVR_CONF_WASAPI_STR  = 'WASAPI';
        BIT_8_INT_DEPT_STR    = '8-bit integer';
        BIT_16_INT_DEPT_STR   = '16-bit integer';
        BIT_24_INT_DEPT_STR   = '24-bit integer';
        BIT_32_FLOAT_DEPT_STR = '32-bit float';
    end
	properties (GetAccess = public, SetAccess = protected)
        dev_name_str     = '';
        divr_conf_str    = '';
        n_chan_val       = -1;
        sam_frq_val      = -1;
        sam_per_fram_val = -1;
        bit_dept_str     = '';
        mic_map_idx_cvec = [];
        cap_tim          = -1;
        max_sam_val      = -1;
        n_fram           = -1;
        wav_mat          = [];
        cap_aud_obj      = NaN;
        sum_n_over_run   = -1;
    end
    methods (Access = public)
        function this = wav_cap_def (new_dev_name_str, ...
                                     new_divr_conf_str, ...
                                     new_n_chan_val, ...
                                     new_sam_frq_val, ...
                                     new_sam_per_fram_val, ...
                                     new_bit_dept_str, ...
                                     new_mic_map_idx_cvec, ...
                                     new_cap_tim)
             this.dev_name_str     = new_dev_name_str;
             this.divr_conf_str    = new_divr_conf_str;
             this.n_chan_val       = new_n_chan_val;
             this.sam_frq_val      = new_sam_frq_val;
             this.sam_per_fram_val = new_sam_per_fram_val;
             this.bit_dept_str     = new_bit_dept_str;
             this.set_cap_tim(new_cap_tim);
             this.mic_map_idx_cvec = new_mic_map_idx_cvec;
        end
        function set_cap_tim (this, new_cap_tim)
            new_max_sam_val  = fix(new_cap_tim * this.sam_frq_val);
            new_max_sam_val  = new_max_sam_val - mod(new_max_sam_val, this.sam_per_fram_val);
            new_n_fram       = fix(new_max_sam_val / this.sam_per_fram_val);
            if (new_max_sam_val <= 0) || (mod(new_max_sam_val, this.sam_per_fram_val) ~= 0)
                error('Error occurred.')
            end
            this.cap_tim     = new_cap_tim;
            this.max_sam_val = new_max_sam_val;
            this.n_fram      = new_n_fram;
        end
        function init (this)
            this.cap_aud_obj = audioDeviceReader ( ...
                'Device',          this.dev_name_str, ...
                'Driver',          this.divr_conf_str, ...
                'NumChannels',     this.n_chan_val, ...
                'SampleRate',      this.sam_frq_val, ...
                'SamplesPerFrame', this.sam_per_fram_val, ...
                'BitDepth',        this.bit_dept_str);
            setup(this.cap_aud_obj);
        end
        function cap (this)
            new_wav_mat = zeros(this.max_sam_val, this.n_chan_val);
            new_sum_n_over_run = 0;
            for i_fram = 1: 1: this.n_fram
                [tmp_frm_cap_aud_mat, n_over_run] = this.cap_aud_obj();
                new_wav_mat(((i_fram - 1) * this.sam_per_fram_val) + 1: (i_fram * this.sam_per_fram_val), :) = tmp_frm_cap_aud_mat;
                new_sum_n_over_run = new_sum_n_over_run + n_over_run;
            end
            this.sum_n_over_run = new_sum_n_over_run;
            this.wav_mat        = new_wav_mat(:, this.mic_map_idx_cvec);
        end
        function save (this, path_str)
            obj           = this; %#ok<NASGU>
            file_name_str = strcat('REC_', datestr(now,'yyyymmddTHHMMSS'), '_CH_', num2str(this.n_chan_val), '.mat');
            save(fullfile(path_str, file_name_str), 'obj');
        end
        function de_init (this)
            release(this.cap_aud_obj);
        end
        function open_setting (this)
            if strcmp(this.divr_conf_str, this.DIVR_CONF_ASIO_STR)
                asiosettings(this.dev_name_str);
            else
                error('Error occurred.')
            end
        end
    end
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
