
function PSTH(trial, n_id, k, neurons_id, smooth)
    % Creates peri-stimulus time histograms (PSTH) for different
    % neural units across different trials. PSTH is the rate as a
    % Spike Density (Average over Several Runs). 
    % Input:
    %        trial: A structure that contains the data. 
    %         n_id: A vector that contains all the trials we want to plot.
    %            k: A number between 1 and 8 that indicates the reaching angle.
    %   neurons_id: A vector that contains all the neuron units we want to plot.
    %       smooth: A number that specifies the number of neighbors that will be
    %               taken around one point to take the average. This will
    %               smooth the histogram to get a continuous rate estimate. 
    % Output:
    %      It creates an histogram.
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
        axis([0,Size,0,1])
        count = zeros(1,Size);
        for n = n_id
            count = count + [trial(n,k).spikes(i,:),zeros(1,Size-size(trial(n,k).spikes,2))];
        end  
        count_smooth = count;
        if smooth > 0
            for s = 1:1:Size
                I = max(1,-smooth+s):1:min(Size,s+smooth);
                for l = I
                    count_smooth(1,s) = count_smooth(1,s)+count(1,l);
                end
                count_smooth(1,s) = (count_smooth(1,s)-count(1,s))/length(I);
            end
        end        
        count_smooth = count_smooth/(0.001*length(n_id));
        bar(1:Size,count_smooth,1.0)
        ylabel({'Neuron',num2str(i)})
    end
    xlabel('time (ms)')
end
