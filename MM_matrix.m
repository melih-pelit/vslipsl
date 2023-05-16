function MM = MM_matrix(X,param)
% MM matrix for testing SLIP models swing leg
% 2019.01.30

%% variables
x_CoM = X(1);
z_CoM = X(2);
theta = X(3);
r_swFoot = X(4);

dx_CoM = X(5);
dz_CoM = X(6);
dtheta = X(7);
dr_swFoot = X(8);

q = [x_CoM; z_CoM; theta; r_swFoot];
dq = [dx_CoM; dz_CoM; dtheta; dr_swFoot];

%% Parameters
L0 = param(1); % [m]
alpha0 = param(2); % angle at which swing leg touches the ground
k0 = param(3); % nominal leg stiffness [N/m]
m_CoM = param(4); % hip mass [kg]

% new parameters for swing leg dynamics
m_swLeg = param(5); % [kg]
m_swFoot = param(6); % [kg]
I_swLeg = param(7); %
I_swFoot = param(8); %
L_thigh = param(9); % [m]
k_swFoot = param(10); %
k_swLeg = param(11); %
theta0 = param(12); %[rad]
r_swFoot0 = param(13); % [m]

g = param(14); % gravitational acc

%% MM Matrix
MM = zeros(4);

MM(1, 1) = m_CoM + m_swLeg + m_swFoot;
MM(1, 2) = 0;
MM(1, 3) = -(sin(theta)*(L_thigh*m_swLeg + 2*L_thigh*m_swFoot + 2*m_swFoot*r_swFoot))/2;
MM(1, 4) = m_swFoot*cos(theta);

MM(2, 1) = 0;
MM(2, 2) = m_CoM + m_swLeg + m_swFoot;
MM(2, 3) = (cos(theta)*(L_thigh*m_swLeg + 2*L_thigh*m_swFoot + 2*m_swFoot*r_swFoot))/2;
MM(2, 4) = m_swFoot*sin(theta);

MM(3, 1) = -(sin(theta)*(L_thigh*m_swLeg + 2*L_thigh*m_swFoot + 2*m_swFoot*r_swFoot))/2;
MM(3, 2) = (cos(theta)*(L_thigh*m_swLeg + 2*L_thigh*m_swFoot + 2*m_swFoot*r_swFoot))/2;
MM(3, 3) = I_swLeg + I_swFoot + (L_thigh^2*m_swLeg)/4 + L_thigh^2*m_swFoot + m_swFoot*r_swFoot^2 + 2*L_thigh*m_swFoot*r_swFoot;
MM(3, 4) = 0;

MM(4, 1) = m_swFoot*cos(theta);
MM(4, 2) = m_swFoot*sin(theta);
MM(4, 3) = 0;
MM(4, 4) = m_swFoot;

end