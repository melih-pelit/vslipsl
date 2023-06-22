function plot_controller_errors(time, ss_controller_info, ds_controller_info, flag_print, varargin)

errors = [
    ss_controller_info.error_ss.Data, ...
    ds_controller_info.error_ds.Data];

errors_size = size(errors);
no_of_plots = errors_size(2);

if nargin == 6
    lgd_font_size = varargin{2};
else
    lgd_font_size = 14;
end

labels = {'h_1 [m]';'h_2 [m]';'h_3 [m]';'h_4 [m]';'h_5 [m/s]'};

fig = figure();
fig.Position = [100 100 800 800]; % make the figure spawn larger

line_width = 1.5;
for i = 1:no_of_plots
    subplot(no_of_plots,1,i)
    plot(time, errors(:,i), 'LineWidth', line_width)
    ylabel(labels(i,:), 'FontSize', lgd_font_size)

    if flag_print
        set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
    end
end

if flag_print
    if nargin >= 5
        export_fig("figures\fig_controller_errors" + varargin{1}, '-m3')
    else
        export_fig("figures\fig_controller_errors", '-m3')
    end
end


