function [ptc, pBest, pBestPtc, vel] = init_ptc(pNumber, f_reset)
%%%%%%%%%% SEED %%%%%%%%%%%%%%%%%%%%%%%%
% PSO parameters
k_swFoot_sd = 2000; %
k_swLeg_sd = 150; %
r_swFoot0_sd = 0.1; % [m] free length of the spring for foot
theta0_sd = pi/2; %[rad] 90 degrees because of the way coordinate systems were set

% PSO initial conditions
dx_CoM_0_sd = 1.5;
dz_CoM_0_sd = 0.98;
dtheta_0_sd = 0;
dr_swFoot_0_sd = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 1;

% ptc(1) = y_M0
% ptc(2) = theta_0
% ptc(3) = r_0
% ptc(4) = dx_M0
% ptc(5) = dtheta_0
% ptc(6) = dr_0
% ptc(7) = k_swLeg
% ptc(8) = k_swFoot
% ptc(9) = theta_s0 % spring free position
% ptc(10) = r_s0 % spring free position

if f_reset == 0
    load('pBest\pBestMemo.mat')
    bestVals = pBestMemo(:,1);
    [val, idx] = max(bestVals);
    pBest = val; % initialzie best fitness value ever
    ptc(1,:) = pBestMemo(idx,2:7);
    pBestPtc = ptc(1,:);
    start_random = 2; % randomize starting from the second particle
else
    % if the program was reset
%     delete pbest/pBestMemo.mat
    pBest = 0;
    pBestPtc = zeros(1, 10);
    pBestMemo = [double(pBest), pBestPtc, zeros(1, 3)];
    save('pbest\pBestMemo', 'pBestMemo')
%     ptc(1,:) = [1364.58488484179,373.831282199266,1.52487722368992,0.174524504094327,2.31899474058299,1.79945998921960,0.00109190116031610,-0.00788255262445179];
    % from 201906191600pBestMemo
%     ptc(1,:) = [1234.48146452840,359.713795854701,1.53820994847663,0.132010371950202,0.00114689211968498,-0.00831580064803621];
%     ptc(1,:) = [1126.76652539497,654.806662180312,1.69409552727176,0.113322333743264,0.000915342798558247,0.0100000000000000];
%     ptc(2,:) = [563.633164264395,109.301454717119,1.40135954052896,-0.153456711922190,0.939569075377513,2.81223855447187,0.00126318939612375,-0.00656469280085745];
    start_random = 1; % randomize starting from the second particle
end

ptc(1,:) = [0.898043902780393,4.47844632885521,0.0575153669902351,3.12550147239918,-0.192359011140840,1.03431700772973,3000.61089745711,4902.54246852384,4.80243294846631,0.0147512151762426];
ptc_preset_size = size(ptc);
start_random = ptc_preset_size + 1;
% rr = normrnd(0,0.5); % the range of randomness when initializing
for i = start_random:pNumber
    
    % (b-a).*rand(1000,1) + a; where b is the upper and a is the lower
    % bound
    ptc(i,1) = (1 - 0.7)*rand() + 0.7; %
    ptc(i,2) = (2*pi - pi)*rand() + pi; 
    ptc(i,3) = (0.3 - (-0.3))*rand() + (-0.3);
    ptc(i,4) = 5*rand();
    ptc(i,5) = (10 - (-10))*rand() + (-10);
    ptc(i,6) = (5 - (-5))*rand() + (-5);
    ptc(i,7) = 5000*rand();
    ptc(i,8) = 2000*rand();
    ptc(i,9) = (2*pi)*rand();
    ptc(i,10) = (2 - (-2))*rand() + (-2);
end
% creating initial velocities
% vel = 2*rand(size(ptc)).*ptc;
vel = normrnd(0,2)*rand(size(ptc)).*ptc; % gaussion distribution with mean 0 sigma 2
end