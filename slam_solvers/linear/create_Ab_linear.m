% CREATE_AB_LINEAR
% 16-831 Fall 2016 - *Stub* Provided
% Computes the A and b matrices for the 2D linear SLAM problem
%
% Arguments: 
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
%     A       - MxN matrix
%     b       - Mx1 vector
%
function [As, b] = create_Ab_linear(odom, obs, sigma_o, sigma_l)


% Useful Constants
n_poses = size(odom, 1) + 1; % +1 for prior on the first pose
n_landmarks = max(obs(:,2));

n_odom = size(odom, 1);
n_obs  = size(obs, 1);

% Dimensions of state variables and measurements (all 2 in this case)
p_dim = 2;
l_dim = 2;
o_dim = size(odom, 2);
m_dim = size(obs(1, 3:end), 2);

% A matrix is MxN, b is Mx1
N = p_dim*n_poses + l_dim*n_landmarks;
M = o_dim*(n_odom+1) + m_dim*n_obs;     % +1 for prior on the first pose

% Initialize matrices
A = zeros(M, N);
b = zeros(M, 1);

% Add odometry and landmark measurements to A, b - including prior on first
% pose
sigma_o_sq = inv(sqrt(sigma_o));
sigma_l_sq = inv(sqrt(sigma_l));

H_o = eye(2);
H_m = -eye(2);

for i = 1:n_poses
    A(2*i-1:2*i,2*i-1:2*i) = sigma_o_sq*H_o;
    if i~=1
        A(2*i-1:2*i,2*i-3:2*i-2) = -sigma_o_sq*H_o;
    end
    
    if i == 1
        b(i:i+1) = [0;0];
    else
        b(2*i-1:2*i) = sigma_o_sq*[odom(i-1,1);odom(i-1,2)];
    end
end

base = p_dim*n_poses;

for i = 1:n_obs
    p_idx = obs(i,1);
    l_idx = obs(i,2);
    %A(base +2*i-1:base+2*i, 2*p_idx-1:2*p_idx) = sigma_l_sq*H_m;
    A(base +2*i-1:base+2*i, 2*p_idx-1:2*p_idx) = sigma_l_sq*H_m;
    A(base +2*i-1:base+2*i, base+2*l_idx-1:base+2*l_idx) = -sigma_l_sq*H_m;
    b(base +2*i-1:base+2*i) = sigma_l_sq*[obs(i,3); obs(i,4)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Your code goes here %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Make A a sparse matrix 
As = sparse(A);