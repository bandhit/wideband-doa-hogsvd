%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name        : Export to PNG File
%
% Description : Export to PNG File
%
% Author      : Bandhit Suksiri
%               Information Systems Engineering
%               Kochi University of Technology
%
% Contact     : bandhit.tni@gmail.com
%
% Logs        : Created: 3 May 2018, Bandhit Suksiri,
%               Updated: 3 May 2018, Bandhit Suksiri.
%
% Copyright 2017,
% Signal Processing & New Generation Network Laboratory (FUKULAB),
% Kochi University of Technology (KUT).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_png (file_name, fig_obj, x_val, Y_val, dpi_str)
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
    
    print(file_name, strcat('-r', dpi_str), '-dpng');
end

% EOF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%