function U_ds = VSLIPSL_DS_controller(x_CoM, X_slipsl_ds, flag, ref_star_ds, slipslParams, gains_VSLIPSL)

y_CoM = X_slipsl_ds(2);
dx_CoM = X_slipsl_ds(3);
dy_CoM = X_slipsl_ds(4);

foot_prev = flag(3) - flag(2);

% SLIP-SL constant params
m_M = slipslParams(1);
m_swLeg = slipslParams(2);
m_swFoot = slipslParams(3);
L_thigh = slipslParams(4);
I_swLeg = slipslParams(5);
I_swFoot = slipslParams(6);
gravi = slipslParams(7);

% SLIP-SL collocation parameters
k0_ss = slipslParams(8);
L0_ss = slipslParams(9);
k0_ds = slipslParams(10);
L0_ds = slipslParams(11);
k_swFoot = slipslParams(12);
k_swLeg = slipslParams(13);
theta0 = slipslParams(14);
r0 = slipslParams(15);

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

ref_output_ds = [y_CoM_star; dx_CoM_star];

A_ds = [
    ((y_CoM - del_y_CoM_star*x_CoM)*(L0_ds - (x_CoM^2 + y_CoM^2)^(1/2)))/((x_CoM^2 + y_CoM^2)^(1/2)*(m_M + m_swLeg + m_swFoot)), (y_CoM*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot) + (del_y_CoM_star*(foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot);
    (x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot),                                                                                                        -((foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot)];

K_ds = [
    del2_y_CoM_star*dx_CoM^2 - K_p_ds*(y_CoM - y_CoM_star) - K_d_ds*(dy_CoM - del_y_CoM_star*dx_CoM) - (k0_ds*y_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1) - gravi*(m_M + m_swLeg + m_swFoot) + k0_ds*y_CoM*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot) - (del_y_CoM_star*(k0_ds*(foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1) - k0_ds*x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1)))/(m_M + m_swLeg + m_swFoot);
    del_dx_CoM_star*dx_CoM - K_v_ds*(dx_CoM - dx_CoM_star) + (k0_ds*(foot_prev - x_CoM)*(L0_ds/(y_CoM^2 + (foot_prev - x_CoM)^2)^(1/2) - 1) - k0_ds*x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot)];

det_A_ds = det(A_ds);
if abs(det_A_ds) >= 4e-9
    U_ds = (A_ds^-1)*K_ds;
else
    U_ds = zeros(2,1);
end

end