function [spikedens] = PSTH(trial, t_start, t_end, n_units, n_trials, direction)

    for jj = direction     
        for j = 1:n_units
            for i = 1:n_trials
                spikecount{j}(i,:) = trial(i,direction(jj)).spikes(j,t_start:t_end-1);
            end
            spikedens{j}(jj,:) = sum(spikecount{j},1)/n_trials;
        end
        
        %uncomment to plot
%         figure()
%  
%         for j = 1:n_units
%             subplot(n_units,1,j)
%             bar(t_start:t_end-1,spikedens(j,:))
%             hold on
%             axis([t_start t_end 0 max(max(spikedens))])
%             xlabel('Time (ms)')
%             ylabel(['Neuron = ' num2str(j)])
%         end

        clear spikecount
    end
    
end