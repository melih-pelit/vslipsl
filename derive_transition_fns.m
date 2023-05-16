% 2020.11.02 Code to find generalized coordianates after the switch from
% SLIP to SLIP-SL in the lift-off (DS to SS)

clear

syms x_M y_M theta r real
syms dx_M dy_M dtheta dr real
syms ddx_M ddy_M ddtheta ddr real
syms m_swLeg m_swFoot real
syms I_swLeg I_swFoot real                                   
syms L_thigh real
syms k0 L0 k_swFoot k_swLeg theta0 r0 real 
syms gravi real 

syms m_M foot_prev real
syms x_CoM y_CoM real
syms dx_CoM dy_CoM real

x_swLeg = x_M + 0.5*L_thigh*cos(theta);
y_swLeg = y_M + 0.5*L_thigh*sin(theta);

x_swFoot = x_M + (L_thigh + r)*cos(theta);
y_swFoot = y_M + (L_thigh + r)*sin(theta);

m_tot = m_M + m_swLeg + m_swFoot;

eqn1 = x_CoM - (m_M*x_M + m_swLeg*x_swLeg + m_swFoot*x_swFoot)/m_tot;
eqn2 = y_CoM - (m_M*y_M + m_swLeg*y_swLeg + m_swFoot*y_swFoot)/m_tot;
eqn3 = x_M + (L_thigh + r)*cos(theta) - foot_prev;
eqn4 = y_M + (L_thigh + r)*sin(theta);

%%
q = [x_M; y_M; theta; r];
dq = [dx_M; dy_M; dtheta; dr];

dx_swLeg = jacobian(x_swLeg,q)*dq;
dy_swLeg = jacobian(y_swLeg,q)*dq;

dx_swFoot = jacobian(x_swFoot, q)*dq;
dy_swFoot = jacobian(y_swFoot, q)*dq;

eqn5 = dx_CoM - (m_M*dx_M + m_swLeg*dx_swLeg + m_swFoot*dx_swFoot)/m_tot
eqn6 = dy_CoM - (m_M*dy_M + m_swLeg*dy_swLeg + m_swFoot*dy_swFoot)/m_tot
eqn7 = dx_swFoot
eqn8 = dy_swFoot

%%
constraint1 = pi < theta;
constraint2 = theta < 2*pi;

solX = solve([eqn1, eqn2, eqn3, eqn4, constraint1, constraint2], [x_M, y_M, theta, r], 'ReturnConditions', true ,'Real',true);
% solX = solve([eqn1, eqn2, eqn3, eqn4], [x_M, y_M, theta, r], 'ReturnConditions', true);

% S = solve([eqn1, eqn2], [x_M, theta],'ReturnConditions',true)

% eqn1 = 

% eqn3 = subs(eqn1, x_M, foot_prev - (L_thigh + r)*cos(theta));
% [soltheta, params, conditions] = solve([eqn3, constraint1, constraint2], theta, 'ReturnConditions',true);
% 
% solx = solve([eqn1, eqn2, constraint1, constraint2], [x_M, theta], 'ReturnConditions', true);
% 
% eqn4 = y_CoM - (m_M*y_M + m_swLeg*y_swLeg + m_swFoot*y_swFoot)/(m_M + m_swLeg + m_swFoot) == 0;
% eqn5 = y_M + (L_thigh + r)*sin(theta) == 0;
% 
% soly = solve([eqn4, eqn5, constraint1, constraint2], [y_M, theta], 'ReturnConditions', true);


%% solutions for x_M, y_M and theta

% calc_x_M = (2*L_thigh*m_M*x_CoM - 2*L_thigh*foot_prev*m_swFoot - L_thigh*foot_prev*m_swLeg + 2*L_thigh*m_swLeg*x_CoM + 2*L_thigh*m_swFoot*x_CoM - 2*foot_prev*m_swFoot*r + 2*m_M*r*x_CoM + 2*m_swLeg*r*x_CoM + 2*m_swFoot*r*x_CoM)/(2*L_thigh*m_M + L_thigh*m_swLeg + 2*m_M*r + 2*m_swLeg*r);
% calc_y_M = (2*y_CoM*(L_thigh + r)*(m_M + m_swLeg + m_swFoot))/(2*L_thigh*m_M + L_thigh*m_swLeg + 2*m_M*r + 2*m_swLeg*r);
% calc_theta = 2*pi - acos((2*(foot_prev - x_CoM)*(m_M + m_swLeg + m_swFoot))/(2*L_thigh*m_M + L_thigh*m_swLeg + 2*m_M*r + 2*m_swLeg*r));

%% calculating velocities after lift-off
% some constants were set (r = 0.3 [m], dr = 0 [m/s])

eqn11 = dx_CoM == dx_M + (dr*m_swFoot*cos(theta))/(m_M + m_swLeg + m_swFoot) - (dtheta*sin(theta)*((L_thigh*m_swLeg)/2 + L_thigh*m_swFoot + m_swFoot*r))/(m_M + m_swLeg + m_swFoot);
eqn22 = dy_CoM == dy_M + (dr*m_swFoot*sin(theta))/(m_M + m_swLeg + m_swFoot) + (dtheta*cos(theta)*((L_thigh*m_swLeg)/2 + L_thigh*m_swFoot + m_swFoot*r))/(m_M + m_swLeg + m_swFoot);
% eqn33 is the time derivative of the constraint (x_m + (L_thigh + r)*cos(theta) = foot_prev)
eqn33 = dx_M + dr*cos(theta) - (L_thigh + r)*sin(theta)*dtheta == 0;

sol_vel = solve([eqn11, eqn22, eqn33], [dx_M, dy_M, dtheta])




