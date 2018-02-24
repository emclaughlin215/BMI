function [] = RasterFunct(trial, numangles, units)

    figure(1)
    for j = 1:numangles
        subplot(numangles, 1, j)
        for i = 1:units
            t = 1:1:length(trial(i,j).spikes(1,:));
            bar(t, trial(i,j).spikes(1,:))
            axis([0 800 0 1])
            xlabel('Time')
            ylabel('Spike')
            hold on
            clear t
        end
    end
    
end