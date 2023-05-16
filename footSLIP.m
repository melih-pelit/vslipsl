function [a, a0, SLIP_flags] = footSLIP(q1, q2, SLIP_flags, SLIP_params, a, a0)
% foot locations of SLIP model to draw on top of 5link model animation

%--- SLIP parameters------------
L0 = SLIP_params(1); % [m]
alpha0 = SLIP_params(2); % angle at which swing leg touches the ground
k0 = SLIP_params(3); % nominal leg stiffness [N/m]
m = SLIP_params(4); % hip mass [kg]
g = SLIP_params(5);
%------------------------------

%--- SLIP flags -----------------
f_SS = SLIP_flags(1);
f_SS_change = SLIP_flags(2);
f_DS_change = SLIP_flags(3);
%------------------------------

%--- calculating spring forces --------
L1 = sqrt((q1 - a0)^2 + q2^2);
L2 = sqrt((q1 - a)^2 + q2^2);

Fs1 = k0 * ((L0/L1) - 1) * [q1 - a0; q2];
Fs2 = k0 * ((L0/L2 )- 1) * [q1 - a; q2];
%------------------------------

if f_SS == 1 && q2 < L0 * sin(alpha0) && f_DS_change == 0 && Fs1(1) >= 0
    f_SS = 0;
    f_SS_change = 0;
    f_DS_change = 1;
    
    a = q1 + L0 * cos(alpha0); % calculating the new foots location
    
    L2 = sqrt((q1 - a)^2 + q2^2);
    Fs2 = k0 * ((L0/L2 )- 1) * [q1 - a; q2];
        
%     display('Double Stance')
elseif f_SS == 0 && Fs1(1) < 0 && f_SS_change == 0
    f_SS = 1;
    f_SS_change = 1;
    f_DS_change = 0;
    a0 = a; % a0 is the previous location of the swing leg touch point
    
    L1 = sqrt((q1 - a0)^2 + q2^2);
    Fs1 = k0 * ((L0/L1) - 1) * [q1 - a0; q2];
%     display('Single Stance')
end

SLIP_flags = [f_SS, f_SS_change, f_DS_change];

end