%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Export to PDF File
%
% Description : Export to PDF File
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 14 September 2017, Bandhit Suksiri,
%               Updated: 14 September 2017, Bandhit Suksiri.
%
% Copyright 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_pdf (file_name, fig_obj, x_val, Y_val)
    x_mgin_val = 0;
    y_mgin_val = 0;
    x_size_val = x_val - 2 * x_mgin_val;
    y_size_val = Y_val - 2 * y_mgin_val;
    set(fig_obj, 'Menubar',             'none');
    set(fig_obj, 'Units','centimeters', 'Position', [0, 0, x_size_val, y_size_val] / 2)
    set(fig_obj, 'PaperUnits',          'centimeters');
    set(fig_obj, 'PaperSize',           [x_val, Y_val]);
    set(fig_obj, 'PaperPosition',       [0, 0, x_size_val, y_size_val]);
    set(fig_obj, 'PaperOrientation',    'portrait');
    
    print(file_name, '-r0', '-dpdf');
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%