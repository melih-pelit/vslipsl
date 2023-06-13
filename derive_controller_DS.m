%% Deriving the DS Controller
clc 
clear

syms x_CoM y_CoM dx_CoM dy_CoM
syms foot foot_prev k0_ds L0_ds
syms m_M m_swLeg m_swFoot gravi

% Equation of motion

L1 = sqrt((x_CoM)^2 + y_CoM^2); % length of the new stance leg
L2 = sqrt((x_CoM - foot_prev)^2 + y_CoM^2); % length of the old stance leg

Fs1 = k0_ds * ((L0_ds/L1) - 1) * [x_CoM; y_CoM]; % Spring force of the new stance leg
Fs2 = k0_ds * ((L0_ds/L2 )- 1) * [x_CoM - foot_prev; y_CoM]; % Spring force of the old stance leg
F = Fs1 + Fs2; % spring force matrix

dq = [dx_CoM; dy_CoM];
ddq = [
    F(1)/(m_M + m_swLeg + m_swFoot);
    (F(2) - (m_M + m_swLeg + m_swFoot)*gravi)/(m_M + m_swLeg + m_swFoot)];

syms u_4 u_5

Fs1_var = (k0_ds + u_4) * ((L0_ds/L1) - 1) * [x_CoM; y_CoM]; % Spring force of the new stance leg
Fs2_var = (k0_ds + u_5) * ((L0_ds/L2 )- 1) * [x_CoM - foot_prev; y_CoM]; % Spring force of the old stance leg
F_var = Fs1_var + Fs2_var; % spring force matrix

ddq_var = [
    F_var(1)/(m_M + m_swLeg + m_swFoot);
    (F_var(2) - (m_M + m_swLeg + m_swFoot)*gravi)/(m_M + m_swLeg + m_swFoot)];

f_ds = [dq; ddq];
% obtaining the u terms
U_ds = [u_4; u_5];
g_ds = [zeros(2,2); simplify(jacobian((ddq_var - ddq), U_ds))]

% test
simplify(f_ds + g_ds*U_ds - [dq; ddq_var]) % =0 means derivation is correct

g_4 = g_ds(:,1);
g_5 = g_ds(:,2);


%% Exact linearization via feedback
syms y_CoM_ref(x_CoM) dx_CoM_ref(x_CoM)
syms K_d_ds K_p_ds K_v_ds

z_ds = [x_CoM; y_CoM; dx_CoM; dy_CoM];

% outputs
h4 = y_CoM - y_CoM_ref;
h5 = dx_CoM - dx_CoM_ref;

% Calculating lie derivatives for DS

Lie_f_ds_h4 = jacobian(h4, z_ds)*f_ds;
Lie2_f_ds_h4 = jacobian(Lie_f_ds_h4, z_ds)*f_ds;

Lie_f_ds_h5 = jacobian(h5,z_ds)*f_ds;

Lie_g4_lie_f_ds_h4 = jacobian(Lie_f_ds_h4, z_ds)*g_4;
Lie_g5_lie_f_ds_h4 = jacobian(Lie_f_ds_h4, z_ds)*g_5;

Lie_g4_h5 = jacobian(h5,z_ds)*g_4;
Lie_g5_h5 = jacobian(h5, z_ds)*g_5;

%% Write to controller_ds.txt
syms y_CoM_star dx_CoM_star
syms del_y_CoM_star del_dx_CoM_star
syms del2_y_CoM_star

fid = fopen('controller_ds.txt', 'wt');

fprintf(fid, 'Derived from derive_controller_DS.m \n');
fprintf(fid, '\n');

% h4
h4 = subs(h4, y_CoM_ref(x_CoM), y_CoM_star);
fprintf(fid, 'h4 = %s;\n', char(h4));

Lie_f_ds_h4 = subs(Lie_f_ds_h4, diff(y_CoM_ref(x_CoM), x_CoM), del_y_CoM_star);
Lie2_f_ds_h4 = subs(Lie2_f_ds_h4, diff(y_CoM_ref(x_CoM), x_CoM, x_CoM), del2_y_CoM_star);
Lie2_f_ds_h4 = subs(Lie2_f_ds_h4, diff(y_CoM_ref(x_CoM), x_CoM), del_y_CoM_star);

fprintf(fid, 'Lie_f_ds_h4 = %s;\n', char(Lie_f_ds_h4));
fprintf(fid, 'Lie2_f_ds_h4 = %s;\n', char(Lie2_f_ds_h4));

Lie_g4_lie_f_ds_h4 = subs(Lie_g4_lie_f_ds_h4, diff(y_CoM_ref(x_CoM), x_CoM), del_y_CoM_star);
Lie_g5_lie_f_ds_h4 = subs(Lie_g5_lie_f_ds_h4, diff(y_CoM_ref(x_CoM), x_CoM), del_y_CoM_star);

fprintf(fid, 'Lie_g4_lie_f_ds_h4 = %s;\n', char(Lie_g4_lie_f_ds_h4));
fprintf(fid, 'Lie_g5_lie_f_ds_h4 = %s;\n', char(Lie_g5_lie_f_ds_h4));

fprintf(fid, '\n');

% h5
h5 = subs(h5, dx_CoM_ref(x_CoM), dx_CoM_star);
Lie_f_ds_h5 = subs(Lie_f_ds_h5, diff(dx_CoM_ref(x_CoM), x_CoM), del_dx_CoM_star);
Lie_g4_h5
Lie_g5_h5

fprintf(fid, 'h5 = %s;\n', char(h5));
fprintf(fid, 'Lie_f_ds_h5 = %s;\n', char(Lie_f_ds_h5));
fprintf(fid, 'Lie_g4_h5 = %s;\n', char(Lie_g4_h5));
fprintf(fid, 'Lie_g5_h5 = %s;\n', char(Lie_g5_h5));

%%
(y_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1))/(m_M + m_swLeg + m_swFoot) - (x_CoM*(L0_ds/(x_CoM^2 + y_CoM^2)^(1/2) - 1)*diff(y_CoM_ref(x_CoM), x_CoM))/(m_M + m_swLeg + m_swFoot)
 


%%



A_ds = [
    Lie_g4_lie_f_ds_h4, Lie_g5_lie_f_ds_h4;
    Lie_g4_h5, Lie_g5_h5];

A_ds = simplify(A_ds)

K_ds = [
    -Lie2_f_ds_h4 - K_d_ds*Lie_f_ds_h4 - K_p_ds*h4;
    -Lie_f_ds_h5 - K_v_ds*h5];

K_ds = simplify(K_ds)


%% Substituting partial derv. /w vars.

A_ds_var = subs(A_ds, diff(y_CoM_ref(x_CoM), x_CoM), del_y_CoM_star)

K_ds_var = subs(K_ds, diff(y_CoM_ref(x_CoM), x_CoM, x_CoM), del2_y_CoM_star);
K_ds_var = subs(K_ds_var, diff(y_CoM_ref(x_CoM), x_CoM), del_y_CoM_star);
K_ds_var = subs(K_ds_var, y_CoM_ref(x_CoM), y_CoM_star);

K_ds_var = subs(K_ds_var, diff(dx_CoM_ref(x_CoM), x_CoM), del_dx_CoM_star);
K_ds_var = subs(K_ds_var, dx_CoM_ref(x_CoM), dx_CoM_star)
