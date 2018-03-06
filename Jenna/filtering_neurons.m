
function out = filtering_neurons(trial, baseline)
% It filters the neurons that we want to use for our analysis.
% We decide if we take into account neuron or not based on if their firing
% rate is higher compared to the baseline firing rate (before movement).
% Input:
%       trial:
%    Baseline: Matrix that contains the baseline firing rate for neurons (row) across all trials for different angles (columns)
% Output:
%       dvcds: 
load('monkeydata_training.mat');

t_movement = 300; % ms
total_time = 400; % ms
for neuron = 1:size(baseline,1)
    for orientation = 1:size(baseline,2)
        t_start_movement = ceil(size(trial(neuron,orientation).spikes,2)/total_time*t_movement + 1);
        spikes_movement = trial(neuron,orientation).spikes(:,t_start_movement:end); % across all trials
        rate_movement_orients(neuron, orientation) = mean(mean(spikes_movement)); % average firing rate for the movement phase for each neuron across all trials
        var_neuron_orient(orientation) = mean(var(spikes_movement')); % across trials for each orientation
    end
    
    % Average across all orientations
    var_movement_neuron(neuron) = mean(var_neuron_orient);
    avg_baseline_neuron(neuron) = mean(baseline(neuron,:)); % 0 - 300ms
    avg_movement_neuron(neuron) = mean(rate_movement_orients(neuron,:)); % 300 - 400ms
    
    % Different measures
    FF(neuron) = var_movement_neuron(neuron)/avg_baseline_neuron(neuron);
    firing_rate_ratio(neuron) = avg_movement_neuron(neuron)/avg_baseline_neuron(neuron);
end

% Decide if we want to get rid of the neuron
for neuron = 1:size(baseline,1)
    if %%FF < 3 % we get rid of neuron ---- %firing_rate_ratio < 3
        trial(neuron,:) = [];
    end
end

end
