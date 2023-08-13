function plot_var_stiff_inputs(time, inputs, flag, param, flag_print, varargin)

%% ----------- parameters ----------------------
L0_ss = param(1); % [m]
k0_ss = param(2); % nominal leg stiffness [N/m]
m_M = param(3); % hip mass [kg]
m_swLeg = param(4); % [kg]
m_swFoot = param(5); % [kg]
I_swLeg = param(6); %
I_swFoot = param(7); %
L_thigh = param(8); % [m]
k_swFoot = param(9); %
k_swLeg = param(10); %
theta0 = param(11); %[rad]
r0 = param(12); % [m]
gravi = param(13); % gravitational acc

L0_ds = param(14); % [m] free length of the springs in the DS
k0_ds = param(15); % [N/m] stiffness of the springs in the DS

bound_cst = param(16); % clipped = max(min(x, upper), lower);

% find the flag change indexes
flag_prev = 1;
state_change_idx = [];
for i = 1:length(time)
    if flag_prev ~= flag(i,1)
        state_change_idx(end+1) = i;
        flag_prev = flag(i,1);
    end
end

if nargin == 7
    lgd_font_size = varargin{2};
else
    lgd_font_size = 14;
end

t_start = 0;
t_end = time(end);

fig = figure();
fig.Position = [100 100 1280 720]; % make the figure spawn larger
nominal = [k0_ss; k_swLeg; k_swFoot; k0_ds; k0_ds];
titles = ["k_{stLeg} + u_1"; "k_{swLeg} + u_2"; "k_{swFoot} + u_3"; "k_{DS} + u_4"; "k_{DS} + u_5"];
y_labels= ["[N/m]"; "[Nm/rad]"; "[N/m]"; "[N/m]"; "[N/m]"];
line_width = 1.5;
for i = 1:5
    subplot(5,1,i)
    plot(time, nominal(i) + inputs(:,i), 'LineWidth', line_width)
    hold on
    plot(time, nominal(i)*ones(length(time),1), 'LineWidth', line_width)
    grid on
    % vline(time(state_change_idx),'r')
    xlim([t_start, t_end])
    title(titles(i,:), 'FontSize', lgd_font_size)
    ylabel(y_labels(i,:), 'FontSize', lgd_font_size)

    if i == 1
        legend('nominal', 'var. stiffness', 'FontSize', lgd_font_size)
    end

    if flag_print
        set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
    end

end

if flag_print
    if nargin >= 6
        export_fig("figures\fig_plot_var_stiff_inputs" + varargin{1}, '-m3')
    else
        export_fig("figures\fig_plot_var_stiff_inputs", '-m3')
    end
end
%%