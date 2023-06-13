% 2021.03.09 -  Creating the SLIP-SL True CoM simulation to add robustness
% via variable stiffness springs.

% Melih Pelit

% clc
clear

%% loading slipsl trajectory (openocl solution)
% load('slipsl_traj\dc_comp2021-03-10-16-20.mat') % of SLIPSLOpenOCL2021-03-09-20-33
% load('slipsl_traj\dc_comp2021-03-18-20-51.mat') % of SLIPSLOpenOCL2021-03-18-20-48
% load('slipsl_traj\dc_comp2021-03-28-16-43.mat')
% load('slipsl_traj\dc_comp2021-03-28-21-07.mat') % traj with cost fns disabled 
% load('slipsl_traj\dc_comp2021-03-29-16-16') % with dy_swFoot(TD) < -0.001 but cost functions disabled
% load('slipsl_traj\dc_comp2021-03-29-22-10') % % with dy_swFoot(TD) < -0.001 of SLIPSLOpenOCL2021-03-29-22-08
% load('slipsl_traj\dc_comp2021-04-01-22-47') % SLIPSLOpenOCL2021-04-01-22-46
% load('slipsl_traj\dc_comp2021-04-28-16-38') % SLIPSLOpenOCL2021-04-28-16-37 !GOOD RESULTS

% load('slipsl_traj\dc_comp2021-06-29-15-15') % SLIPSLOpenOCL2021-06-29-15-04
% load('slipsl_traj\dc_comp2021-06-30-15-27') % SLIPSLOpenOCL2021-06-30-15-22

% load('slipsl_traj\dc_comp2021-07-08-17-31') % SLIPSLOpenOCL2021-07-08-17-31

% 2021.09.01: new SLIP-SL trajectory with low CoT (CoT = 0.1907)
load('slipsl_traj\dc_comp2021-09-01-16-06') % SLIPSLOpenOCL2021-09-01-15-36 (@C:\Matlab Workspace\SLIPSL_OpenOCL\202108111747 - SLIPSL_OpenOCL (True CoM))

% optimized SLIP-SL parameters
L0_ss = dc.col_param.L0_ss;
L0_ds = dc.col_param.L0_ds;

k0_ss = dc.col_param.k0_ss;
k0_ds = dc.col_param.k0_ds;

k_swFoot = dc.col_param.k_swFoot;
k_swLeg = dc.col_param.k_swLeg;

theta_s0 = dc.col_param.theta0;
r_s0 = dc.col_param.r0;

dc.simout_ss = dc.simout_ss';
dc.simout_ds = dc.simout_ds';

%% initializing

% constant parameters
% L0 = 1; % [m]
% alpha0 = 1.2; % angle at which swing leg touches the ground
m_M = dc.const_param.m_M; % hip mass [kg]
m_swLeg = dc.const_param.m_swLeg; % [kg]
m_swFoot = dc.const_param.m_swFoot; % [kg]
L_thigh = dc.const_param.L_thigh; % [m]
I_swLeg = dc.const_param.I_swLeg; %
I_swFoot = dc.const_param.I_swFoot; %
g = dc.const_param.gravi; % gravitational acc

% initial condinition
x_M_init = dc.simout_ss(1,1);
y_M_init = dc.simout_ss(1,2);
theta_init = dc.simout_ss(1,3);
r_init = dc.simout_ss(1,4);

x_swLeg_init = x_M_init + (L_thigh/2)*cos(theta_init);
x_swFoot_init = x_M_init + (L_thigh + r_init)*cos(theta_init);
x_CoM_ss_init = (m_M*x_M_init + m_swLeg*x_swLeg_init + m_swFoot*x_swFoot_init)/(m_M + m_swLeg + m_swFoot);

dx_M_init = dc.simout_ss(1,5);
dy_M_init = dc.simout_ss(1,6);
dtheta_init = dc.simout_ss(1,7);
dr_init = dc.simout_ss(1,8);

% initialzating flags
%%%%% initial foot position %%%%%
% foot0 = mean(liftOff.foot - liftOff.x);
% foot0 = 0; % it starts from VLO

% flag = [f_ss; foot; foot_prev; x_M_LO]
flag_ref_ob = 0;
init_flag = [1; dc.col_param.foot; -dc.col_param.footPlus; x_M_init; x_CoM_ss_init; flag_ref_ob]; % start with single stance phase (init_flag(1) == 1)
init_f_record = 0;
param = [L0_ss; k0_ss; m_M; m_swLeg; m_swFoot; I_swLeg; I_swFoot; L_thigh; k_swFoot; k_swLeg; theta_s0; r_s0; g; L0_ds; k0_ds];

x_CoM_init = (x_M_init*m_M + (x_M_init + (L_thigh/2)*cos(theta_init))*m_swLeg + (x_M_init + (L_thigh + r_init)*cos(theta_init))*m_swFoot)/(m_M + m_swLeg + m_swFoot);
y_CoM_init = (y_M_init*m_M + (y_M_init + (L_thigh/2)*sin(theta_init))*m_swLeg + (y_M_init + (L_thigh + r_init)*sin(theta_init))*m_swFoot)/(m_M + m_swLeg + m_swFoot);

dx_CoM_init = (dx_M_init*m_M + (dx_M_init - (L_thigh/2)*sin(theta_init)*dtheta_init)*m_swLeg ...
    + (dx_M_init + dr_init*cos(theta_init) - (L_thigh + r_init)*sin(theta_init)*dtheta_init)*m_swFoot)/(m_M + m_swLeg + m_swFoot);
dy_CoM_init = (dy_M_init*m_M + (dy_M_init + (L_thigh/2)*cos(theta_init)*dtheta_init)*m_swLeg ...
    + (dy_M_init + dr_init*sin(theta_init) + (L_thigh + r_init)*cos(theta_init)*dtheta_init)*m_swFoot)/(m_M + m_swLeg + m_swFoot);

init_cond = [x_M_init; y_M_init; theta_init; r_init; x_CoM_init; y_CoM_init; dx_M_init; dy_M_init; dtheta_init; dr_init; dx_CoM_init; dy_CoM_init];


%% Loading the referece trajectory

% load('ref_traj_star\ref_star2021-04-23-15-31.mat') % dc_comp2021-04-01-22-47
% load('ref_traj_star\ref_star2021-04-27-21-28.mat') % dc_comp2021-03-29-16-16
% load('ref_traj_star\ref_star2021-04-28-16-43') % dc_comp2021-04-28-16-38 !! GOOD RESULT

% load('ref_traj_star\ref_star2021-06-29-15-30') % dc_comp2021-06-29-15-15
% load('ref_traj_star\ref_star2021-06-30-15-28') % dc_comp2021-06-30-15-27
% load('ref_traj_star\ref_star2021-06-30-16-52') % dc_comp2021-06-29-15-15, CoM based references

% load('ref_traj_star\ref_star2021-07-09-14-56') % dc_comp2021-07-08-17-31, increasing dy_sw(TD) constraint to <= -0.1

% 2021.09.01: new SLIP-SL trajectory with low CoT (CoT = 0.1907)
load('ref_traj_star\ref_star2021-09-03-14-20') % dc_comp2021-09-01-16-06

ref_star_ss = [ref_star.ss.x_M_star, ... 
    ref_star.ss.y_M_star, ref_star.ss.x_swFoot_star, ref_star.ss.y_swFoot_star, ...
    ref_star.ss.del_y_M_star, ref_star.ss.del_x_swFoot_star, ref_star.ss.del_y_swFoot_star, ...
    ref_star.ss.del2_y_M_star, ref_star.ss.del2_x_swFoot_star, ref_star.ss.del2_y_swFoot_star, ...
    ref_star.ss.x_CoM_star];

ref_star_ds = [ref_star.ds.x_CoM_star, ...
    ref_star.ds.y_CoM_star, ref_star.ds.dx_CoM_star, ...
    ref_star.ds.del_y_CoM_star, ref_star.ds.del_dx_CoM_star, ...
    ref_star.ds.del2_y_CoM_star];


%% Controller gains
% Controller Gains

gain.K_p = 200;
gain.K_d = 40;

gain.K_p_sw = 200;
gain.K_d_sw = 40;

gain.K_p_ds = 200;
gain.K_d_ds = 40;
gain.K_v_ds = 40; % 5 in visser 2012

% gain.K_p = 40;
% gain.K_d = 15;
% 
% gain.K_p_sw = 40;
% gain.K_d_sw = 15;
% 
% gain.K_p_ds = 40;
% gain.K_d_ds = 15;
% gain.K_v_ds = 15; % 5 in visser 2012

gains = [gain.K_p; gain.K_d; gain.K_p_sw; gain.K_d_sw; gain.K_p_ds; gain.K_d_ds; gain.K_v_ds];
%%
% Run Simulation
flag_dist = 0; % disturbance flag

Tf = 10; % final time
sample_time = 0.001;

% open_system('slip_w_sw_leg_DEtest')
% sim('slip_w_sw_leg_DEtest')

open_system('vslipsl_sim')
sim('vslipsl_sim')

%%
f_record = 0;
if f_record == 1
    filename = sprintf('pBestMemo%s.mat', datestr(now,'yyyy-mm-dd-HH-MM'));
    subfolder = 'pBest';
    save(fullfile(subfolder,filename),'pBestMemo')
end

%% Plot
% find the flag change indexes
flag_prev = 1;
state_change_idx = [];
for i = 1:length(time)
    if flag_prev ~= flag(i,1)
        state_change_idx(end+1) = i;
        flag_prev = flag(i,1);
    end
end

t_start = 0;
t_end = time(end);

figure()
nominal = [k0_ss; k_swLeg; k_swFoot; k0_ds; k0_ds];
y_labels = ["u_1"; "u_2"; "u_3"; "u_4"; "u_5"];
for i = 1:5
    subplot(5,1,i)
    plot(time, nominal(i) + inputs(:,i))
    hold on
    plot(time, nominal(i)*ones(length(time),1))
    grid on
    % vline(time(state_change_idx),'r')
    xlim([t_start, t_end])
    ylabel(y_labels(i,:))

    if i == 4
        ylim([16000, 18500]);
    end
end

% debug
figure()
subplot(5,1,1)
plot(time, ss_controller_info.det_A_ss.Data)
ylabel('det(A_{ss})')
subplot(5,1,2)
plot(time, ss_controller_info.K_ss.Data)
ylabel('K_{ss}')
subplot(5,1,3)
plot(time, ss_controller_info.error_ss.Data)
ylabel('error_{ss}')
subplot(5,1,4)
plot(time, ss_controller_info.Lie_f_h.Data)
ylabel('Lie_f_h')
subplot(5,1,5)
plot(time, ss_controller_info.Lie2_f_h.Data)
ylabel('Lie2_f_h')
%% Animation

f_animation = 0; % animation flag
f_video = 0; % to turn on video recording
frame_leap = 10;
f_pause = 0;

if f_animation == 1
    animation(simout, flag, param, time, frame_leap, f_video, f_pause, sample_time)
end

%% Analysis
f_analysis = 0; % flag for turning on and of the analysis
if f_analysis == 1
    
    x_CoM = simout(:, 1);
    z_CoM = simout(:, 2);
    dz_CoM  = simout(:, 6);
    
    % z_CoM plot
    figure
    plot(time, z_CoM)
    xlabel('time [s]')
    ylabel('z_CoM [m]')
    
    % Limit cycle
    figure
    plot(z_CoM, dz_CoM)
    
    % Simple Poincare Map (z_CoM vs vel ratio)
    figure
    hold on
    plot(qVLO(:,2), qVLO(:,6)./qVLO(:,5))
    plot(qVLO(:,2), qVLO(:,6)./qVLO(:,5), '.b')
    xlabel('$z_{CoM}$','Interpreter','Latex');
    ylabel('$\dot{z}_{CoM}/\dot{x}_{CoM}$','Interpreter','Latex')
    
    % Mean distance between VLO points
    figure
    diff = qVLO(1:(length(qVLO)-1), :) - qVLO(2:length(qVLO), :);
%     for i = 1:length(diff)
%         dist(i) = mean(diff(i,2:8));
%     end
%     plot(1:length(diff), dist)
    
    % Finding the stationary point]
%     L = length(qVLO);
%     last = round(L*0.1);
%     fixed = [mean(qVLO(last:L,1)),mean(qVLO(last:L,2)),mean(qVLO(last:L,3)),mean(qVLO(last:L,4)),mean(qVLO(last:L,5)),mean(qVLO(last:L,6)),mean(qVLO(last:L,7)),mean(qVLO(last:L,8))];
%     
    % checking the trajectories
    
    %%%%% finding beginning of the ss phase %%%%%%%
    i = length(time) - round(length(time)*0.01); % Start Search
    if flag(i,1) == 1
        % if randomly chosen location is single stance
        while (1)
            if flag(i,1) == 0
                break
            end
            i = i + 1;
        end
        while (1)
            if flag(i,1) == 1
                break
            end
            i = i + 1;
        end
    else
        % if randomly chosen location is double stance
        while (1)
            if flag(i,1) == 1
                break
            end
            i = i + 1;
        end
    end
    startP = i;
    while (1)
        if flag(i,1) == 0
            break
        end
        i = i + 1;
    end
    while (1)
        if flag(i,1) == 1
            break
        end
        i = i + 1;
    end
    endP = i - 1; % end of the double stance phase
    %%%%%%%%%%%
    
    %%%%% assigning the variables %%%%%
    slipsl.x = simout(startP:endP,1);
    slipsl.z = simout(startP:endP,2);
    slipsl.theta = simout(startP:endP,3);
    slipsl.r_swFoot = simout(startP:endP,4);
    slipsl.dx = simout(startP:endP,5);
    slipsl.dz = simout(startP:endP,6);
    slipsl.dtheta = simout(startP:endP,7);
    slipsl.dr_swFoot = simout(startP:endP,8);
    slipsl.footX = slipsl.x + (L_thigh + slipsl.r_swFoot).*cos(2*pi - slipsl.theta);
    slipsl.footZ = slipsl.z + (L_thigh + slipsl.r_swFoot).*sin(2*pi - slipsl.theta);
    slipsl.footDZ = slipsl.dz + slipsl.dr_swFoot.*sin(2*pi - slipsl.theta) + (L_thigh + slipsl.r_swFoot).*cos(2*pi - slipsl.theta).*slipsl.dtheta;
    slipsl.footFoot = flag(startP:endP,2);
    slipsl.time = time(startP:endP);
    %%%%%%%%%%
    
    figure
    subplot(3,2,1);
    plot(slipsl.time, slipsl.x, 'LineWidth', 2)
    grid on
    title('x comparison')
    
    subplot(3,2,2);
    plot(slipsl.time, slipsl.dx, 'LineWidth', 2)
    grid on
    title('dx comparison')
    
    subplot(3,2,3);
    plot(slipsl.time, slipsl.z, 'LineWidth', 2)
    grid on
    title('z comparison')
    
    subplot(3,2,4);
    plot(slipsl.time, slipsl.dz, 'LineWidth', 2)
    grid on
    title('dz comparison')
    
    subplot(3,2,5);
    plot(slipsl.time, slipsl.theta, 'LineWidth', 2)
    grid on
    title('theta comparison')
    
    subplot(3,2,6);
    plot(slipsl.time, slipsl.r_swFoot, 'LineWidth', 2)
    grid on
    title('r_{swFoot}')
    
    figure
    subplot(2,2,1)
    plot(slipsl.time, slipsl.footX, 'LineWidth', 2)
    grid on;
    title('foot X location comparison')
    
    subplot(2,2,2)
    plot(slipsl.time, slipsl.footZ, 'LineWidth', 2)
    grid on
    title('foot Z location comparison')
    
    subplot(2,2,3:4)
    plot(slipsl.time, slipsl.footDZ, 'LineWidth', 2)
    grid on
    title('foot dz  comparison')
end

