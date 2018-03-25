% MEAS_LANDMARK_JACOBIAN
% 16-831 Fall 2016 - *Stub* Provided
% Compute the Jacobian of the measurement function
%
% Arguments: 
%     rx    - robot's x position
%     ry    - robot's y position
%     lx    - landmark's x position
%     ly    - landmark's y position
%
% Returns:
%     H     - Jacobian of the measurement fuction
%
function H = meas_landmark_jacobian(rx, ry, lx, ly)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Your code goes here %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = (lx-rx)^2+(ly-ry)^2;
H_m_p = [(ly-ry)/t, (rx-lx)/t;(rx-lx)/sqrt(t), (ry-ly)/sqrt(t)];
H_m_l = [(ry-ly)/t, (lx-rx)/t;(lx-rx)/sqrt(t), (ly-ry)/sqrt(t)];
H = [H_m_p, H_m_l];