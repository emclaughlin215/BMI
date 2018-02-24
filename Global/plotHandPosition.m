
function plotHandPosition(trial, n_id, k)
    % Creates a figure for hand positions for different trials.
    % Input:
    %   trial: A structure that contains the data. 
    %    n_id: A vector that contains all the trials we want to plot.
    %       k: A number between 1 and 8 that indicates the reaching angle.
    % Output:
    %      It creates the hand position plot.
    figure
    subplot(3,1,1) % x
    hold on
    for n = n_id
        plot(trial(n, k).handPos(1,:))
    end
    title('x direction')
    
    subplot(3,1,2) % y
    hold on
    for n = n_id
        plot(trial(n, k).handPos(2,:))
    end
    title('y direction')
    
    subplot(3,1,3) % z
    hold on
    for n = n_id
        plot(trial(n, k).handPos(3,:))
    end
    title('z direction')
end
