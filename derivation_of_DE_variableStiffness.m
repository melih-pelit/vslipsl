% 2021/04/06 - Derivation of the dynamic equations for SLIP-SL with
% variable stiffness springs (Single Stance Phase)

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

%%
% velocities
dx_swLeg = jacobian(x_swLeg, q)*dq;
dy_swLeg = jacobian(y_swLeg, q)*dq;

dx_swFoot = jacobian(x_swFoot, q)*dq;
dy_swFoot = jacobian(y_swFoot, q)*dq;


%%

KE_M = 0.5*m_M*(dx_M^2 + dy_M^2); % kinetic energ of the main mass
KE_swLeg = 0.5*m_swLeg*(dx_swLeg^2 + dy_swLeg^2) + 0.5*I_swLeg*dtheta^2;
KE_swFoot = 0.5*m_swFoot*(dx_swFoot^2 + dy_swFoot^2) + 0.5*I_swFoot*dtheta^2;

KE = KE_M + KE_swLeg + KE_swFoot;
PE = m_M*gravi*y_M + m_swLeg*gravi*y_swLeg + m_swFoot*gravi*y_swFoot...
    + 0.5*(k0_ss)*(L0_ss - sqrt((x_M - foot)^2 + y_M^2))^2 + 0.5*(k_swLeg)*(theta0 - theta)^2 + 0.5*(k_swFoot)*(r0 - r)^2;

% Potential energy of SLIP-SL with variable stifnesses
PE_var = m_M*gravi*y_M + m_swLeg*gravi*y_swLeg + m_swFoot*gravi*y_swFoot...
    + 0.5*(k0_ss + u1)*(L0_ss - sqrt((x_M - foot)^2 + y_M^2))^2 + 0.5*(k_swLeg + u2)*(theta0 - theta)^2 + 0.5*(k_swFoot + u3)*(r0 - r)^2;

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

%%
% dz = [dq; calc_ddq_var];
u = [u1; u2; u3];

u_terms = simplify(tau_var - tau);
S = simplify(jacobian(u_terms, u)); % adding 4 zeros because dq is not a effected directly by u_i

test1 = simplify(tau_var - M*ddq - C - G -S*u) % this being 0 confirms that (tau_var = M*ddq + C + G +g*u)

f = [dq; simplify(inv(M)*(- C - G))];
g = [zeros(4,3); simplify(-inv(M)*S)];

% test2 = simplify([dq; inv(M)*( - C - G - S*u)] - (f + g*u))

%%
g1 = g(:,1);
g2 = g(:,2);
g3 = g(:,3);

%% Relative Degree (Isidori pg. 235 (i) )
z = [q; dq];

syms y_M_ref(x_M)

h1 = y_M - y_M_ref;

Lie_g1_h1 = jacobian(h1, z)*g1; % =0
Lie_g2_h1 = jacobian(h1, z)*g2; % =0
Lie_g3_h1 = jacobian(h1, z)*g3; % =0

Lie_f_h1 = jacobian(h1, z)*f;

Lie_g1_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g1); % != 0
Lie_g2_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g2); % != 0
Lie_g3_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g3); % != 0

% rel. degree of h1 is 2

%%
syms theta_ref(x_M)

h2 = theta - theta_ref;

Lie_g1_h2 = jacobian(h2, z)*g1; % =0
Lie_g2_h2 = jacobian(h2, z)*g2; % =0
Lie_g3_h2 = jacobian(h2, z)*g3; % =0

Lie_f_h1 = jacobian(h1, z)*f;

Lie_g1_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g1); % != 0
Lie_g2_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g2); % != 0
Lie_g3_lie_f_h1 = simplify(jacobian(Lie_f_h1, z)*g3); % != 0

%% CoM Position

x_CoM = simplify((m_M*x_M + m_swLeg*x_swLeg + m_swFoot*x_swFoot)/(m_M + m_swLeg + m_swFoot));
y_CoM = simplify((m_M*y_M + m_swLeg*y_swLeg + m_swFoot*y_swFoot)/(m_M + m_swLeg + m_swFoot));

dx_CoM = simplify(jacobian(x_CoM, q)*dq);
dy_CoM = simplify(jacobian(y_CoM, q)*dq);

ddx_CoM = simplify(jacobian(dx_CoM, q_dq)*dq_ddq);
ddy_CoM = simplify(jacobian(dy_CoM, q_dq)*dq_ddq);