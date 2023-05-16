% 2020/10/28 - Derivation of dynamic equation for SLIP model with passive
% swing leg x-y coordianate notation is used

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

q = [x_M; y_M; theta; r];
dq = [dx_M; dy_M; dtheta;  dr];
ddq = [ddx_M; ddy_M; ddtheta;  ddr];

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
    + 0.5*k0_ss*(L0_ss - sqrt((x_M - foot)^2 + y_M^2))^2 + 0.5*k_swLeg*(theta0 - theta)^2 + 0.5*k_swFoot*(r0 - r)^2;

%%
% % writing kinetic energy from OneNote notes
% KE1 = 0.5*m_swLeg*((dx_CoM + L_thigh*0.5*dtheta*cos(theta + pi/2))^2 + (dz_CoM + L_thigh*0.5*dtheta*sin(theta + pi/2))^2) ... 
%     + 0.5*I_swLeg*dtheta^2 ...
%     + 0.5*m_swFoot*((dx_CoM + dr_swFoot*cos(theta) + (L_thigh + r_swFoot)*dtheta*cos(theta + pi/2))^2 + (dz_CoM + dr_swFoot*sin(theta) + (L_thigh + r_swFoot)*dtheta*sin(theta + pi/2))^2) ...
%     + 0.5*I_swFoot*dtheta^2;
% 
% KE = KE1 + 0.5*m_CoM*(dx_CoM^2 + dz_CoM^2); % adding the KE of CoM
% 
% % potential energy
% PE1 = m_swLeg*gravi*(z_CoM - L_thigh*0.5*sin(theta)) + m_swFoot*gravi*(z_CoM - (L_thigh + r_swFoot)*sin(theta)) + 0.5*k_swLeg*(theta0 - theta)^2 + 0.5*k_swFoot*(r_swFoot0 - r_swFoot)^2;
% 
% PE = PE1 + m_CoM*gravi*z_CoM; % adding the PE of CoM

% simplifying the energy functions
KE = simplify(KE);
PE = simplify(PE);
DE = 0;

% calculating the Lagrangian
Lag = simplify(KE - PE);
q_dq = [q; dq];
dq_ddq = [dq; ddq];

%% Calculating the dynamical equation
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

calc_ddq = simplify(inv(M)*(-C-G))

%% CoM Position

x_CoM = simplify((m_M*x_M + m_swLeg*x_swLeg + m_swFoot*x_swFoot)/(m_M + m_swLeg + m_swFoot));
y_CoM = simplify((m_M*y_M + m_swLeg*y_swLeg + m_swFoot*y_swFoot)/(m_M + m_swLeg + m_swFoot));

dx_CoM = simplify(jacobian(x_CoM, q)*dq);
dy_CoM = simplify(jacobian(y_CoM, q)*dq);

ddx_CoM = simplify(jacobian(dx_CoM, q_dq)*dq_ddq);
ddy_CoM = simplify(jacobian(dy_CoM, q_dq)*dq_ddq);

%%
% EoM: M*ddq + C + G = 0 --> ddq = inv(M)*(-C-G)

% % pring the functions
% fprintf('-----------inertia_matrix-------------\n\n')
% for i=1:1:length(q)
%     for j= 1:1:length(q)
%         fprintf('M%d%d', i ,j)
%         M(i,j)
%     end
% end
% fprintf('-----------Coliori and centrifugal force matrix (vector)-------------\n\n')
% for i=1:1:length(q)
%     fprintf('C%d', i)
%     C(i,:)
% end
%  fprintf('-----------gravitational_vector-----------\n\n')
% for i=1:1:length(q)
%     fprintf('G%d', i)
%     G(i,:)
% end
% save('calculation_result.txt','M','C','G','KE_all','PE_all')