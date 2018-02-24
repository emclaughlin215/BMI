function [] = PSTH(trial, t_start, t_end, n_units, n_trials, direction)
    %t = 1:(t_end-t_start);
    for jj = direction
        
        for j = 1:n_units
            for i = 1:n_trials
                spikecount{j}(i,:) = trial(i,direction(jj)).spikes(j,t_start:t_end-1);
            end
            spikedens(j,:) = sum(spikecount{j},1)/n_trials;
        end

        maxYSpike = max(max(spikedens));
        
        figure()
 
        for j = 1:n_units
            subplot(n_units,1,j)
            bar(t_start:t_end-1,spikedens(j,:))
            hold on
            axis([t_start t_end 0 maxYSpike])
            xlabel('Time (ms)')
            ylabel(['Neuron = ' num2str(j)])
        end
        clear spikecount
        clear spikedens
    end
    
end