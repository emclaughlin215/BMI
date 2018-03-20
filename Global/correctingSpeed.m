
function newSpeed = correctingSpeed(Param, x, v)
% It corrects the velocity by taking into account the current velocity
% and the preferred direction by using a control system.
% Inputs:
%   Param: Parameters of the network.
%   x: A vector with the current position.
%   v: A vector with the curent velocity.
% Output:
%   newSpeed: A vector with the corrected velocity.

    preferred_direction = [cos(Param.prefdir), sin(Param.prefdir)];
%     preferred_direction = [cos(degtorad(120)), sin(degtorad(120))];
%     x = [-0.10 0.05];
%     v = [-0.3 0.3];
%         
    % Find the corrected speed using a control system
    k_v = 0.5; % control gain: importance to difference in velocity
    k_x = 0.01; % control gain: importance to difference in position
          
    diff_position = pdist([x; preferred_direction], 'euclidean'); % Stronger effect the further away you are
    error = preferred_direction - v; % between desired and actual velocity
    newSpeed = k_v * error + k_x * diff_position;
    
    figure
    hold on
    plot(x(1),x(2),'rx')
    quiver(0, 0, preferred_direction(1), preferred_direction(2),'m')
    quiver(x(1), x(2), v(1), v(2),'b')
    quiver(x(1), x(2), newSpeed(1), newSpeed(2),'g')
    legend('Current position', 'Preferred direction', 'Current velocity', 'New velocity')
end