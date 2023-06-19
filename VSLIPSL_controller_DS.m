function [U_ds, error_ds, det_A_ds, reference_traj] = VSLIPSL_controller_DS(x_CoM, X_slipsl_ds, flag, ref_star_ds, param, gains_VSLIPSL)

y_CoM = X_slipsl_ds(2);
dx_CoM = X_slipsl_ds(3);
dy_CoM = X_slipsl_ds(4);

foot_prev = flag(3) - flag(2);

% SLIP-SL constant params
% SLIP-SL constant params
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

% Controller Gains
K_p_ds = gains_VSLIPSL(5);
K_d_ds = gains_VSLIPSL(6);
K_v_ds = gains_VSLIPSL(7); % 5 in visser 2012

% interpolation
x_CoM_star_ds = ref_star_ds(:,1);
y_CoM_star = interp1(x_CoM_star_ds, ref_star_ds(:,2), x_CoM);
dx_CoM_star = interp1(x_CoM_star_ds, ref_star_ds(:,3), x_CoM);

del_y_CoM_star = interp1(x_CoM_star_ds, ref_star_ds(:,4), x_CoM);
del_dx_CoM_star = interp1(x_CoM_star_ds, ref_star_ds(:,5), x_CoM);

del2_y_CoM_star = interp1(x_CoM_star_ds, ref_star_ds(:,6), x_CoM);

% below are derived using derive_controller_DS.m
% which outputs to ds_controller.txt
h4 = y_CoM - y_CoM_star;
Lie_f_ds_h4 = dy_CoM - del_y_CoM_star*dx_CoM;
Lie2_f_ds_h4 = (k0_ds*y_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1) - gravi*(m_M + m_swLeg + m_swFoot) + k0_ds*y_CoM*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot) - del2_y_CoM_star*dx_CoM^2 + (del_y_CoM_star*(k0_ds*(foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1) - k0_ds*x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1)))/(m_M + m_swLeg + m_swFoot);
Lie_g4_lie_f_ds_h4 = (y_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot) - (del_y_CoM_star*x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot);
Lie_g5_lie_f_ds_h4 = (y_CoM*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot) + (del_y_CoM_star*(foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot);

h5 = dx_CoM - dx_CoM_star;
Lie_f_ds_h5 = - del_dx_CoM_star*dx_CoM - (k0_ds*(foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1) - k0_ds*x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot);
Lie_g4_h5 = (x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot);
Lie_g5_h5 = -((foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot);

A_ds = [
    Lie_g4_lie_f_ds_h4, Lie_g5_lie_f_ds_h4;
    Lie_g4_h5, Lie_g5_h5];

K_ds = [
    -Lie2_f_ds_h4 - K_d_ds*Lie_f_ds_h4 - K_p_ds*h4;
    -Lie_f_ds_h5 - K_v_ds*h5];

det_A_ds = det(A_ds);
if abs(det_A_ds) >= 4e-9
    U_ds = (A_ds^-1)*K_ds;
else
    U_ds = zeros(2,1);
end

reference_traj = [y_CoM_star; dx_CoM_star];
error_ds = [h4; h5];

end