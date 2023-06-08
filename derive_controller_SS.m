%% Derivation of controller for V-SLIP-SL with h2 and h3 as the
% swing foot x and y position

clear all
close all

syms x_M y_M theta r real
syms dx_M dy_M dtheta dr real
syms ddx_M ddy_M ddtheta ddr real
syms m_swLeg m_swFoot real
syms I_swLeg I_swFoot real                                   
syms L_thigh real
syms k0_ss L0_ss k_swFoot k_swLeg theta0 r0 real 
syms gravi real 

syms m_M foot

syms u1 u2 u3

q = [x_M; y_M; theta; r];
dq = [dx_M; dy_M; dtheta;  dr];
ddq = [ddx_M; ddy_M; ddtheta; ddr];

% positions of M, swLeg and sw foot

% x position of M is x_M and y positon of M is y_M

x_swLeg = x_M + 0.5*L_thigh*cos(theta);
y_swLeg = y_M + 0.5*L_thigh*sin(theta);

x_swFoot = x_M + (L_thigh + r)*cos(theta);
y_swFoot = y_M + (L_thigh + r)*sin(theta);

%% velocities
dx_swLeg = jacobian(x_swLeg, q)*dq;
dy_swLeg = jacobian(y_swLeg, q)*dq;

dx_swFoot = jacobian(x_swFoot, q)*dq;
dy_swFoot = jacobian(y_swFoot, q)*dq;

KE_M = 0.5*m_M*(dx_M^2 + dy_M^2); % kinetic energ of the main mass
KE_swLeg = 0.5*m_swLeg*(dx_swLeg^2 + dy_swLeg^2) + 0.5*I_swLeg*dtheta^2;
KE_swFoot = 0.5*m_swFoot*(dx_swFoot^2 + dy_swFoot^2) + 0.5*I_swFoot*dtheta^2;

KE = KE_M + KE_swLeg + KE_swFoot;
PE = m_M*gravi*y_M + m_swLeg*gravi*y_swLeg + m_swFoot*gravi*y_swFoot...
    + 0.5*(k0_ss)*(L0_ss - sqrt(x_M^2 + y_M^2))^2 + 0.5*(k_swLeg)*(theta0 - theta)^2 + 0.5*(k_swFoot)*(r0 - r)^2;

% Potential energy of SLIP-SL with variable stifnesses
PE_var = m_M*gravi*y_M + m_swLeg*gravi*y_swLeg + m_swFoot*gravi*y_swFoot...
    + 0.5*(k0_ss + u1)*(L0_ss - sqrt(x_M^2 + y_M^2))^2 + 0.5*(k_swLeg + u2)*(theta0 - theta)^2 + 0.5*(k_swFoot + u3)*(r0 - r)^2;

%% Calculating the dynamical equation
% simplifying the energy functions
KE = simplify(KE);
PE = simplify(PE);
DE = 0;

% calculating the Lagrangian
Lag = simplify(KE - PE);
q_dq = [q; dq];
dq_ddq = [dq; ddq];

for i = 1:1:length(q)
Ldq(i,:) = diff(Lag,dq(i,:));
Ldq_dt(i,:) = jacobian(Ldq(i,:),q_dq)*dq_ddq;
Lq(i,:) = diff(Lag,q(i,:));
DEdq(i,:) = diff(DE,dq(i,:));
tau(i,:)= Ldq_dt(i,:) - Lq(i,:) + DEdq(i,:);
end

% Getting the matrices of DE
M = simplify(jacobian(tau,ddq));
G = simplify(collect(tau-subs(tau,gravi,0),gravi));
M_line = M*ddq;
C = simplify(tau - M_line - G);

H = simplify(tau - M_line);

calc_ddq = simplify(inv(M)*(-C-G));

%% Calculating the dynamical equation (Variable stiffness)
% simplifying the energy functions
KE = simplify(KE);
PE_var = simplify(PE_var);
DE = 0;

% calculating the Lagrangian
Lag = simplify(KE - PE_var);
q_dq = [q; dq];
dq_ddq = [dq; ddq];

for i = 1:1:length(q)
Ldq(i,:) = diff(Lag,dq(i,:));
Ldq_dt(i,:) = jacobian(Ldq(i,:),q_dq)*dq_ddq;
Lq(i,:) = diff(Lag,q(i,:));
DEdq(i,:) = diff(DE,dq(i,:));
tau_var(i,:)= Ldq_dt(i,:) - Lq(i,:) + DEdq(i,:);
end

%% dz = [dq; calc_ddq_var];
u = [u1; u2; u3];

u_terms = simplify(tau_var - tau);
S = simplify(jacobian(u_terms, u)); % adding 4 zeros because dq is not a effected directly by u_i

test1 = simplify(tau_var - M*ddq - C - G -S*u) % this being 0 confirms that (tau_var = M*ddq + C + G +g*u)

f = [dq; simplify(inv(M)*(- C - G))];
g = [zeros(4,3); simplify(-inv(M)*S)];

% test2 = simplify([dq; inv(M)*( - C - G - S*u)] - (f + g*u))

%% CoM Position
x_CoM = simplify((m_M*x_M + m_swLeg*x_swLeg + m_swFoot*x_swFoot)/(m_M + m_swLeg + m_swFoot));
y_CoM = simplify((m_M*y_M + m_swLeg*y_swLeg + m_swFoot*y_swFoot)/(m_M + m_swLeg + m_swFoot));

dx_CoM = simplify(jacobian(x_CoM, q)*dq);
dy_CoM = simplify(jacobian(y_CoM, q)*dq);

ddx_CoM = simplify(jacobian(dx_CoM, q_dq)*dq_ddq);
ddy_CoM = simplify(jacobian(dy_CoM, q_dq)*dq_ddq);

% sw foot position and velocity

x_sw = x_M + (L_thigh + r)*cos(theta);
y_sw = y_M + (L_thigh + r)*sin(theta);

dx_sw = jacobian(x_sw, q)*dq
dy_sw = jacobian(y_sw, q)*dq

g1 = g(:,1);
g2 = g(:,2);
g3 = g(:,3);

%% Relative Degree (Isidori pg. 235 (i) )
z = [q; dq];

syms y_M_ref(x_M)

h1 = y_M - y_M_ref;

Lie_g1_h1 = simplify(jacobian(h1, z)*g1); % =0
Lie_g2_h1 = simplify(jacobian(h1, z)*g2); % =0
Lie_g3_h1 = simplify(jacobian(h1, z)*g3); % =0

Lie_f_h1 = simplify(jacobian(h1, z)*f);

Lie_g1_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g1); % != 0
Lie_g2_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g2); % != 0
Lie_g3_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g3); % != 0

% rel. degree of h1 is 2

syms x_swFoot_ref(x_M)

x_swFoot = x_M + (L_thigh + r)*cos(theta);
h2 = x_swFoot - x_swFoot_ref;

Lie_g1_h2 = simplify(jacobian(h2, z)*g1) % =0
Lie_g2_h2 = simplify(jacobian(h2, z)*g2) % =0
Lie_g3_h2 = simplify(jacobian(h2, z)*g3) % =0

Lie_f_h2 = simplify(jacobian(h2, z)*f);

Lie_g1_lie_f_h2 = simplify(jacobian(Lie_f_h2, z)*g1) % != 0
Lie_g2_lie_f_h2 = simplify(jacobian(Lie_f_h2, z)*g2) % != 0
Lie_g3_lie_f_h2 = simplify(jacobian(Lie_f_h2, z)*g3) % != 0

% relative degree of h2 is 2


%% h3 Relative Degree
syms y_swFoot_ref(x_M)

y_swFoot = y_M + (L_thigh + r)*sin(theta);
h3 = y_swFoot - y_swFoot_ref;

Lie_g1_h3 = simplify(jacobian(h3, z)*g1) % =0
Lie_g2_h3 = simplify(jacobian(h3, z)*g2) % =0
Lie_g3_h3 = simplify(jacobian(h3, z)*g3) % =0

Lie_f_h3 = simplify(jacobian(h3, z)*f);

Lie_g1_lie_f_h3 = simplify(jacobian(Lie_f_h3, z)*g1) % != 0
Lie_g2_lie_f_h3 = simplify(jacobian(Lie_f_h3, z)*g2) % != 0
Lie_g3_lie_f_h3 = simplify(jacobian(Lie_f_h3, z)*g3) % != 0

Lie2_f_h1 = simplify(jacobian(Lie_f_h1,z)*f);
Lie2_f_h2 = simplify(jacobian(Lie_f_h2,z)*f);
Lie2_f_h3 = simplify(jacobian(Lie_f_h3,z)*f);


%% Substituting the partial derivatives with variables
syms y_M_star x_swFoot_star y_swFoot_star
syms del_y_M_star del_x_swFoot_star del_y_swFoot_star % these are all partial derivatives wrt. x_M
syms del2_y_M_star del2_x_swFoot_star del2_y_swFoot_star % these are all second order partial derivatives wrt. x_M

fid = fopen('controller_ss.txt', 'wt');

fprintf(fid, 'Derived from derive_controller_SS.m \n');
fprintf(fid, '\n');

% h1
h1 = subs(h1, y_M_ref(x_M), y_M_star)

Lie_g1_h1 
Lie_g2_h1 
Lie_g3_h1

Lie_f_h1 = subs(Lie_f_h1, diff(y_M_ref(x_M), x_M), del_y_M_star)

Lie_g1_lie_f_h1 = subs(Lie_g1_lie_f_h1, diff(y_M_ref(x_M), x_M), del_y_M_star)
Lie_g2_lie_f_h1 = subs(Lie_g2_lie_f_h1, diff(y_M_ref(x_M), x_M), del_y_M_star)
Lie_g3_lie_f_h1 = subs(Lie_g3_lie_f_h1, diff(y_M_ref(x_M), x_M), del_y_M_star)

fprintf(fid, 'h1 = %s;\n', char(h1));

fprintf(fid, 'Lie_g1_h1 = %s;\n', char(Lie_g1_h1));
fprintf(fid, 'Lie_g2_h1 = %s;\n', char(Lie_g2_h1));
fprintf(fid, 'Lie_g3_h1 = %s;\n', char(Lie_g3_h1));

fprintf(fid, 'Lie_f_h1 = %s;\n', char(Lie_f_h1));

fprintf(fid, 'Lie_g1_lie_f_h1 = %s;\n', char(Lie_g1_lie_f_h1));
fprintf(fid, 'Lie_g2_lie_f_h1 = %s;\n', char(Lie_g2_lie_f_h1));
fprintf(fid, 'Lie_g3_lie_f_h1 = %s;\n', char(Lie_g3_lie_f_h1));

fprintf(fid, '\n');

% h2

h2 = subs(h2, x_swFoot_ref(x_M), x_swFoot_star)

Lie_g1_h2
Lie_g2_h2
Lie_g3_h2

Lie_f_h2 = subs(Lie_f_h2, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star)

Lie_g1_lie_f_h2 = subs(Lie_g1_lie_f_h2, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star)
Lie_g2_lie_f_h2 = subs(Lie_g2_lie_f_h2, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star)
Lie_g3_lie_f_h2 = subs(Lie_g3_lie_f_h2, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star)

fprintf(fid, 'h2 = %s;\n', char(h2));

fprintf(fid, 'Lie_g1_h2 = %s;\n', char(Lie_g1_h2));
fprintf(fid, 'Lie_g2_h2 = %s;\n', char(Lie_g2_h2));
fprintf(fid, 'Lie_g3_h2 = %s;\n', char(Lie_g3_h2));

fprintf(fid, 'Lie_f_h2 = %s;\n', char(Lie_f_h2));

fprintf(fid, 'Lie_g1_lie_f_h2 = %s;\n', char(Lie_g1_lie_f_h2));
fprintf(fid, 'Lie_g2_lie_f_h2 = %s;\n', char(Lie_g2_lie_f_h2));
fprintf(fid, 'Lie_g3_lie_f_h2 = %s;\n', char(Lie_g3_lie_f_h2));

fprintf(fid, '\n');

% h3

h3 = subs(h3, y_swFoot_ref(x_M), y_swFoot_star);

Lie_g1_h3
Lie_g2_h3
Lie_g3_h3

Lie_f_h3 = subs(Lie_f_h3, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);

Lie_g1_lie_f_h3 = subs(Lie_g1_lie_f_h3, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);
Lie_g2_lie_f_h3 = subs(Lie_g2_lie_f_h3, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);
Lie_g3_lie_f_h3 = subs(Lie_g3_lie_f_h3, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);

fprintf(fid, 'h3 = %s;\n', char(h3));

fprintf(fid, 'Lie_g1_h3 = %s;\n', char(Lie_g1_h3));
fprintf(fid, 'Lie_g2_h3 = %s;\n', char(Lie_g2_h3));
fprintf(fid, 'Lie_g3_h3 = %s;\n', char(Lie_g3_h3));

fprintf(fid, 'Lie_f_h3 = %s;\n', char(Lie_f_h3));

fprintf(fid, 'Lie_g1_lie_f_h3 = %s;\n', char(Lie_g1_lie_f_h3));
fprintf(fid, 'Lie_g2_lie_f_h3 = %s;\n', char(Lie_g2_lie_f_h3));
fprintf(fid, 'Lie_g3_lie_f_h3 = %s;\n', char(Lie_g3_lie_f_h3));

fprintf(fid, '\n');

Lie2_f_h1 = subs(Lie2_f_h1, diff(y_M_ref(x_M), x_M, x_M), del2_y_M_star);
Lie2_f_h1 = subs(Lie2_f_h1, diff(y_M_ref(x_M), x_M), del_y_M_star);

Lie2_f_h2 = subs(Lie2_f_h2, diff(x_swFoot_ref(x_M), x_M, x_M), del2_x_swFoot_star);
Lie2_f_h2 = subs(Lie2_f_h2, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star);

Lie2_f_h3 = subs(Lie2_f_h3, diff(y_swFoot_ref(x_M), x_M, x_M), del2_y_swFoot_star);
Lie2_f_h3 = subs(Lie2_f_h3, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);

fprintf(fid, 'Lie2_f_h1 = %s;\n', char(Lie2_f_h1));
fprintf(fid, 'Lie2_f_h2 = %s;\n', char(Lie2_f_h2));
fprintf(fid, 'Lie2_f_h3 = %s;\n', char(Lie2_f_h3));

fclose(fid);


%% Deriving the controller (SS)
% clear y_M_ref theta_ref r_ref

A_ss = [
    Lie_g1_lie_f_h1, Lie_g2_lie_f_h1, Lie_g3_lie_f_h1;
    Lie_g1_lie_f_h2, Lie_g2_lie_f_h2, Lie_g3_lie_f_h2;
    Lie_g1_lie_f_h3, Lie_g2_lie_f_h3, Lie_g3_lie_f_h3];

A_ss = simplify(A_ss)

syms K_d K_p K_d_sw K_p_sw

K_ss = [
    -Lie2_f_h1 - K_d*Lie_f_h1 - K_p*h1;
    -Lie2_f_h2 - K_d*Lie_f_h2 - K_p*h2;
    -Lie2_f_h3 - K_d*Lie_f_h3 - K_p*h3];

K_ss = simplify(K_ss)


%% Substituting the partial derivatives with variables
syms y_M_star x_swFoot_star y_swFoot_star
syms del_y_M_star del_x_swFoot_star del_y_swFoot_star % these are all partial derivatives wrt. x_M
syms del2_y_M_star del2_x_swFoot_star del2_y_swFoot_star % these are all second order partial derivatives wrt. x_M

A_ss_var = subs(A_ss, diff(y_M_ref(x_M), x_M), del_y_M_star);
A_ss_var = subs(A_ss_var, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star);
A_ss_var = subs(A_ss_var, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);

K_ss_var = subs(K_ss, diff(y_M_ref(x_M), x_M, x_M), del2_y_M_star);
K_ss_var = subs(K_ss_var, diff(x_swFoot_ref(x_M), x_M, x_M), del2_x_swFoot_star);
K_ss_var = subs(K_ss_var, diff(y_swFoot_ref(x_M), x_M, x_M), del2_y_swFoot_star);

K_ss_var = subs(K_ss_var, diff(y_M_ref(x_M), x_M), del_y_M_star);
K_ss_var = subs(K_ss_var, diff(x_swFoot_ref(x_M), x_M), del_x_swFoot_star);
K_ss_var = subs(K_ss_var, diff(y_swFoot_ref(x_M), x_M), del_y_swFoot_star);

K_ss_var = subs(K_ss_var, y_M_ref(x_M), y_M_star);
K_ss_var = subs(K_ss_var, x_swFoot_ref(x_M), x_swFoot_star);
K_ss_var = subs(K_ss_var, y_swFoot_ref(x_M), y_swFoot_star);

simplify(det(A_ss_var))


%% Adding regularization term to the contoller
syms eps_1 eps_2 eps_3
syms u_1 u_2 u_3
syms inv_A_ss_var % it takes a lot of time to take the actual inverse

eps = diag([eps_1, eps_2, eps_3])
E = ((eps + eye(3,3))^-1)


det_A_ss_var = det(A_ss_var)




