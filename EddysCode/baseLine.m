function [averageSpikeDens] = baseLine(trial)
    
    %parameters for baseline calc
    n_trials = 100;
    n_units = 98;
    direction = 1:8;
    t_start = 1;
    t_end = 300;
    
    for jj = direction
        for j = 1:n_units
            for i = 1:n_trials
                spikecount{j,jj}(i,:) = trial(i,direction(jj)).spikes(j,t_start:t_end-1);
            end
            spikedens{j,jj}(1,:) = sum(spikecount{j,jj},1)/n_trials;
            averageSpikeDens(j,jj) = mean(spikedens{j,jj});
        end
        clear spikecount
        clear spikedens
    end
end