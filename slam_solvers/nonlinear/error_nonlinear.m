% ERROR_NONLINEAR
% 16-831 Fall 2016 - *Stub* Provided
% Computes the total error of all measurements (odometry and landmark)
% given the current state estimate
%
% Arguments: 
%     x       - Current estimate of the state vector
%     odom    - Matrix that contains the odometry measurements
%               between consecutive poses. Each row corresponds to
%               a measurement. 
%                 odom(:,1) - x-value of odometry measurement
%                 odom(:,2) - y-value of odometry measurement
%     obs     - Matrix that contains the landmark measurements and
%               relevant information. Each row corresponds to a
%               measurement.
%                 obs(:,1) - idx of pose at which measurement was 
%                   made
%                 obs(:,2) - idx of landmark being observed
%                 obs(:,3) - x-value of landmark measurement
%                 obs(:,4) - y-value of landmark measurement
%     sigma_o - Covariance matrix corresponding to the odometry
%               measurements
%     sigma_l - Covariance matrix corresponding to the landmark
%               measurements
% Returns:
%     err     - total error of all measurements
%
function err = error_nonlinear(x, odom, obs, sigma_odom, sigma_landmark)
%% Extract useful constants which you may wish to use
n_poses = size(odom, 1) + 1;                % +1 for prior on the first pose
n_landmarks = max(obs(:,2));

n_odom = size(odom, 1);
n_obs  = size(obs, 1);

% Dimensions of state variables and measurements (all 2 in this case)
p_dim = 2;                                  % pose dimension
l_dim = 2;                                  % landmark dimension
o_dim = size(odom, 2);                      % odometry dimension
m_dim = size(obs(1, 3:end), 2);    % landmark measurement dimension

% A matrix is MxN, b is Nx1
N = p_dim*n_poses + l_dim*n_landmarks;
M = o_dim*(n_odom+1) + m_dim*n_obs;         % +1 for prior on the first pose

%% Initialize error
err = 0;
b = zeros(M, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Your code goes here %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_o_sq = inv(sqrt(sigma_odom));
sigma_l_sq = inv(sqrt(sigma_landmark));


for i = 1:n_poses
    if i == 1
        b(i:i+1) = [0;0];
    else
        rx1 = x(2*i-3);
        ry1 = x(2*i-2);
        rx2 = x(2*i-1);
        ry2 = x(2*i);
        b(2*i-1:2*i) = -[odom(i-1,1);odom(i-1,2)] + meas_odom(rx1, ry1, rx2, ry2);
    end
end

base = p_dim*n_poses;

for i = 1:n_obs
    p_idx = obs(i,1);
    l_idx = obs(i,2);

    rx = x(2*p_idx-1);
    ry = x(2*p_idx);
    lx = x(base+2*l_idx-1);
    ly = x(base+2*l_idx);
    
    delta = [obs(i,3); obs(i,4)] - meas_landmark(rx, ry, lx, ly);
    b(base +2*i-1:base+2*i) = -delta;
end

tmp = b;
for i = base+1:2:length(tmp)
    if tmp(i) > pi || tmp(i) < -pi
        tmp(i) = wrapToPi(tmp(i));
    end
end

err = err + tmp'*tmp;
