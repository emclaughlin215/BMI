
function [direction, angle] = getPopulationVector(trial,n,k)
    % Finds the global neuronal population vector for a trial for a reching angle.
    % It is used to predict arm movements.
    % Input:
    %   trial: A structure that contains the data. 
    %       n: A number that is the trial we want to plot.
    %       k: A number between 1 and 8 that indicates the reaching angle.
    % Output:
    %      direction: A vector that contains x and y direction of final movement.
    %          angle: The angle (in degrees) of the final movement.
    I = size(trial(1,1).spikes,1);
    [~, ~, vector] = tuning(trial, 1:I, 0);
    
    Weights = sum(trial(n,k).spikes,2)/size(trial(n,k).spikes,2);
    Weights = Weights/sum(Weights);
    
    direction = Weights'*vector;
    direction_norm = direction/sqrt(direction(1,1)^2+direction(1,2)^2);
    angle = atan(direction(1,2)/direction(1,1))*180/pi;
    
    figure
    plot([0,direction_norm(1,1)],[0,direction_norm(1,2)])
    axis([-1,1,-1,1])
    title('Population vector')
end
