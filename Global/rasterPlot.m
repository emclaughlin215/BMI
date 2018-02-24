
function rasterPlot(trial, n_id, k, neurons_id)
    % Creates raster plot (time on x-axis and neural units in y-axis) for
    % different neural units across different trials.
    % For each neuron unit, we plot the spikes of whatever number of trials we want.
    % Input:
    %        trial: A structure that contains the data. 
    %         n_id: A vector that contains all the trials we want to plot.
    %            k: A number between 1 and 8 that indicates the reaching angle.
    %   neurons_id: A vector that contains all the neuron units we want to plot.
    % Output:
    %      It creates a figure for the raster plot.
    color = rand(length(n_id),3);
    sizes = zeros(length(n_id),1);
    for n = n_id
        sizes(n-n_id(1,1)+1,1) = size(trial(n,k).spikes,2);
    end
    Size = max(sizes);
    
    figure;
    for i = neurons_id
        subplot(length(neurons_id),1,i-neurons_id(1,1)+1)
        if i == neurons_id(1,1)
            title(['Trials ',num2str(n_id(1,1)),' to ',num2str(n_id(1,end))])
        end
        hold on
        axis([0,Size,0,1])
        for n = n_id            
            for j = 1:1:size(trial(n,k).spikes,2)
                if trial(n,k).spikes(i,j) == 1
                    plot([j,j],[0,1],'Color',color(n-n_id(1,1)+1,:))
                end
            end
        end    
        hold off
        ylabel({'Neuron',num2str(i)})
    end
    xlabel('time (ms)')
end
