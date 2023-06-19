function plot_controller_errors(time, ss_controller_info, ds_controller_info, flag_print, varargin)

errors = [
    ss_controller_info.error_ss.Data, ...
    ds_controller_info.error_ds.Data];

errors_size = size(errors);
no_of_plots = errors_size(2);

labels = ['h_1';'h_2';'h_3';'h_4';'h_5'];

figure()
for i = 1:no_of_plots
    subplot(no_of_plots,1,i)
    plot(time, errors(:,i))
    ylabel(labels(i,:))

    if flag_print
        set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
    end
end

if flag_print
    if nargin == 5
        export_fig("figures\fig_controller_errors" + varargin{1}, '-m3')
    else
        export_fig("figures\fig_controller_errors", '-m3')
    end
end


