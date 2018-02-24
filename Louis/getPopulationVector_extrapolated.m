function [direction, angle] = getPopulationVector_extrapolated(trial,n,k)
    I=size(trial(1,1).spikes,1);
    [~, ~, vector, ~] = tuning_extrapolated(trial, 1:I, 0);
    
    Weights = sum(trial(n,k).spikes,2)/size(trial(n,k).spikes,2);
    Weights = Weights/sum(Weights);
    
    direction = Weights'*vector;
    direction_norm = direction/sqrt(direction(1,1)^2+direction(1,2)^2);
    angle = atan(direction(1,2)/direction(1,1))*180/pi;
    plot([0,direction_norm(1,1)],[0,direction_norm(1,2)])
    axis([-1,1,-1,1])
    title('Population vector')
    
end