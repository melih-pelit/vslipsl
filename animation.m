function animation(simout, flag, param, time, frame_leap, f_video, f_pause, sample_time, start_end_time)
figure()
simout_size = size(simout);
nt = simout_size(1);

% ----------- parameters ----------------------
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

% x_CoM calculated


if f_video == 1
    video_v = VideoWriter('SLIP_w_sw_leg.avi');
    open(video_v);
end

for i = (start_end_time(1)/sample_time):frame_leap:(start_end_time(2)/sample_time)
    
    X = simout(i,:);
    
    % variables
    x_M = X(1);
    y_M = X(2);
    theta = X(3);
    r = X(4);
    x_CoM = X(5);
    y_CoM = X(6);
    
    
    dx_M = X(7);
    dy_M = X(8);
    dtheta = X(9);
    dr = X(10);
    dx_CoM = X(11);
    dy_CoM = X(12);
    
    f_SS = flag(i,1);
    foot = flag(i,2);
    foot_prev = flag(i,3);
    
    
    if f_SS == 1
        % Single stance phase
        % locations of mass and the link
        
        hip = [x_M; y_M];
        knee = [x_M + L_thigh*cos(theta); y_M + L_thigh*sin(theta)];
        foot_sw = [x_M + (L_thigh + r)*cos(theta); y_M + (L_thigh + r)*sin(theta)];
        plot(hip(1), hip(2), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
        
        hold on
                
        plot([hip(1), knee(1)], [hip(2), knee(2)], 'r', 'LineWidth', 2);
        plot(foot_sw(1), foot_sw(2), 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
        plot([hip(1), foot], [hip(2), 0], 'g', 'LineWidth', 2) % plotting the stance leg
        
        plot(x_CoM, y_CoM, 'o')
           
    else
        % Double stance phase
        plot(x_CoM, y_CoM, 'o', 'MarkerSize', 7, 'MarkerFaceColor', 'b');
        hold on
        plot([x_CoM, foot], [y_CoM, 0], 'g', 'LineWidth', 2) % plotting the stance leg
        plot([x_CoM, foot_prev], [y_CoM, 0], 'r', 'LineWidth', 2) % plotting the swing leg
        
    end
    
    plot([-1000, 1000], [0, 0], 'k') % plotting the ground
    txt1 = ['f_{SS}: ' num2str(f_SS)];
    text(-1, 2.5, txt1)
    
    xlim([-2 + x_CoM, 2 + x_CoM])
    ylim([-0.5, 1.5])
    title(['t = ' num2str(time(i))])
    hold off
    
    if f_video == 1
        frame = getframe(gcf);
        writeVideo(video_v,frame);
    end
    
    if f_pause == 1
        pause
    else
        pause(0.01)
    end
    
end
end