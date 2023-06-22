function plot_trajectory_tracking(time, simout, flag, param, ss_controller_info, ds_controller_info, flag_print, varargin)

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

%%
% Flag
foot = flag(:,2);

% States
x_M = simout(:,1);
y_M = simout(:,2);
theta = simout(:,3);
r = simout(:,4);
x_CoM = simout(:,5);
y_CoM = simout(:,6);

dx_CoM = simout(:,11);

% 
x_swFoot = x_M + (L_thigh + r).*cos(theta) - foot;
y_swFoot = y_M + (L_thigh + r).*sin(theta);

y_M(flag(:,1) ~= 1) = NaN;
x_swFoot(flag(:,1) ~= 1) = NaN;
y_swFoot(flag(:,1) ~= 1) = NaN;

simout_selected_ss = [y_M, x_swFoot, y_swFoot];
simout_selected_ds = [y_CoM, dx_CoM];

simout_selected = [simout_selected_ss, simout_selected_ds];

ref_traj = [
    ss_controller_info.reference_traj_ss.Data, ...
    ds_controller_info.reference_traj_ds.Data];
ref_traj(1,:) = NaN;
size_ref_traj = size(ref_traj);
num_of_plots = size_ref_traj(2);

if nargin == 9
    lgd_font_size = varargin{2};
else
    lgd_font_size = 14;
end

labels = [
    "$y_M$ [m]"; 
    "$x_{swFoot}$ [m]"; 
    "$y_{swFoot}$ [m]"; 
    "$y_{CoM}$ [m]";
    "$\dot{x}_{CoM}$ [m/s]"];

fig = figure();
fig.Position = [100 100 800 800]; % make the figure spawn larger
line_width = 1.5;
tcl = tiledlayout(5,1);
for i = 1:num_of_plots
%     subplot(num_of_plots,1,i)
    nexttile(tcl)
    plot(time, simout_selected(:, i), 'LineWidth', line_width, 'Color', [0.8500 0.3250 0.0980])
    hold on
    plot(time, ref_traj(:, i), "k--", 'LineWidth', 1.0)
    

    ylabel(labels(i,:), Interpreter="latex", FontSize=lgd_font_size)

    if i == num_of_plots
        xlabel('Time [sec]', 'FontSize', lgd_font_size)
    end

    if i == 1
        legend('refence', 'actual', 'FontSize', lgd_font_size)
    end

    if flag_print
        set(gcf, 'color', 'none');
        set(gca, 'color', 'none');
    end
end

if flag_print
    if nargin >= 8
        export_fig("figures\fig_trajectory_tracking" + varargin{1}, '-m3')
    else
        export_fig("figures\fig_trajectory_tracking", '-m3')
    end
end