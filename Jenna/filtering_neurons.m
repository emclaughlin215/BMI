
function out = filtering_neurons(trial, type)
% It filters the neurons that we want to use for our analysis.
% We decide if we take into account the neuron or not based on if their firing
% rate is higher compared to the baseline firing rate (before movement).
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles
%       type: A string that can be the type of filtering we want to do (FF, Partial)
% Output:
%        out: A structure that contains baseline, rate and speed (divided into bins) for each orientation. 

    neural_data = getNeuronData(trial);
    baseline = neural_data(1,1).baseline;
    
    num_angles = size(trial,2); % number angles
    num_trials = size(trial,1); % number trials

    t_movement = 300; % ms
    total_time = 500; % ms
    for i = 1:size(baseline,1) % neuron
        for k = 1:size(baseline,2) % orientation
            t_start_movement = ceil(size(trial(i,k).spikes,2)/total_time*t_movement + 1);
            spikes_movement = trial(i,k).spikes(:,t_start_movement:end); % across all trials
            rate_movement_orients(i, k) = mean(mean(spikes_movement)); % average firing rate for the movement phase for each neuron across all trials
            var_neuron_orient(k) = mean(var(spikes_movement')); % variance across trials for each orientation
        end

        % Average variance across all orientations and across all trials for each neural unit
        var_movement_neuron(i,1) = mean(var_neuron_orient);
        
        % Average firing rate across all directions and across all trials for each neural unit
        avg_baseline_neuron(i,1) = mean(baseline(i,:)); % 0 - 300ms
        avg_movement_neuron(i,1) = mean(rate_movement_orients(i,:)); % 300 - 500ms

        % Different measures
        FF(i,1) = var_movement_neuron(i)/avg_baseline_neuron(i);
        firing_rate_ratio(i,1) = avg_movement_neuron(i)/avg_baseline_neuron(i);
    end
    
    % Decide if we want to get rid of the neuron
    for neuron = 1:size(baseline,1)
        switch type
            case 'FF'
                if FF(neuron) < 3 % we get rid of neuron
                    bool_neurons(neuron,1) = 0; % we remove
                else
                    bool_neurons(neuron,1) = 1; % we keep
                end
        end
    end
    bool_neurons = logical(bool_neurons);
    
    % Remove rate and baseline of this neuron for all
    % orientations and all trials
    for k = 1:num_angles
        for n = 1:num_trials
            neural_data(n,k).PSTH(~bool_neurons,:) = [];
        end
    end
    avg_baseline_neuron(~bool_neurons,:) = [];

    % Find the longest velocity
    % Initialize
    for k = 1:num_angles
        max_length_vel(k) = -inf;
        for n = 1:num_trials
            if length(neural_data(n,k).handvel(1,:)) > max_length_vel(k)
                max_length_vel(k) = length(neural_data(n,k).handvel(1,:));
            end
        end
    end
            
    % Average speed across all trials for each direction
    for k = 1:num_angles
        v_x = NaN(num_trials,max_length_vel(k));
        v_y = NaN(num_trials,max_length_vel(k));
        for n = 1:num_trials
            current_length = length(neural_data(n,k).handvel(1,:));
            v_x(n,1:current_length) = neural_data(n,k).handvel(1,:);
            v_y(n,1:current_length) = neural_data(n,k).handvel(2,:);             
        end
        avg_velocity{k} = [nanmean(v_x); nanmean(v_y)];
    end
    
    % Average PSTH across all trials for each direction for each neuron
    % All PSTH same for each iteration???
    for k = 1:num_angles
        rate{k} = neural_data(1,k).PSTH;
    end
    
    bins = 5; % number of divisions we want
    [rate_split, velocity_split] = data_stepping(trial, bins, rate, avg_velocity);
 
    % Output
    for k = 1:size(baseline,2)
        out(k).baseline = avg_baseline_neuron;
        out(k).rate = rate_split{k};
        out(k).speed = velocity_split{k};
        out(k).bool_neurons = bool_neurons;
    end
end

function [rate_split, velocity_split] = data_stepping(trial, bins, rate, avg_velocity)
% Splits the average firing rate and velocity (x and y) into the average
% of each bin for every orientation.
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
%       bins: Number of divisions of data.
%       rate: A struct with the PSTH for all neural units across time for each orientation.
%       avg_velocity: A struct with the velocity (x and y) across time for each orientation.
% Output:
%       rate_split: A struct with the PSTH for all neural units across bins for each orientation
%       avg_velocity: A struct with the velocity (x and y) across bins for each orientation.

    num_angles = size(trial,2); % select number of angles to consider
    num_trials = size(trial,1);
    
    % Rate
    for k = 1:num_angles
        T = size(rate{k},2); % time
        bin_count = floor(T/bins);
        for num_bin = 1:bins
            rate_split{k}(:,num_bin) = ...
                sum(rate{k}(:,(num_bin-1)*bin_count+1:num_bin*bin_count),2)/bin_count;
        end
    end
    
    % Velocity
    for k = 1:num_angles
        T = size(avg_velocity{k},2); % time
        bin_count = floor(T/bins);
        for num_bin = 1:bins
            velocity_split{k}(:,num_bin) = ...
                sum(avg_velocity{k}(:,(num_bin-1)*bin_count+1:num_bin*bin_count),2)/bin_count;
        end
    end
end

function neural_data = getNeuronData(trial)
% Calculates hand velocity, hand position, spikes without baseline, 
% PSTH and baseline for every trial in each direction.
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
% Output:
%       neural_data: A structure that contains the PSTH, hand position
%                    and hand velocity for each trial for each direction.
    numangles = 8; % select number of angles to consider
    
    %% Baseline
    % Parameters for baseline function
    params_baseline.n_trials = 100; % number trials
    params_baseline.n_units = 98; % number units
    params_baseline.t_start = 1; % start time
    params_baseline.t_end = 300; % end time
    params_baseline.direction = 1:numangles; % directions of movement

    % Calculate the baseline of each neuron
    baseline_spikedens = baseLine(trial, params_baseline);

    %% PSTH
    % Parameters for PSTH function
    params_PSTH.n_trials = 50; % average over this number of trials
    params_PSTH.n_units = 98; % number units
    params_PSTH.t_start = 300; % start time
    params_PSTH.t_end = 500; % end time
    params_PSTH.direction = 1:numangles; % directions of movement

    % Calculate the PSTH according to parameters chosen above
    spikedens = PSTH(trial, params_PSTH);

    %% Hand velocity
    neural_data = handVelocity(trial, spikedens, baseline_spikedens, numangles, params_baseline);
end

%% Hand velocity function

function neural_data = handVelocity(trial, spikedens, baseline_spikedens, numangles, params)
% This calculates the hand velocity.
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
%       spikedens: For each neural unit, we return the average spike rate
%                  at each ms of movement for each direction.
%       baseline_spikedens: Average spike rate of the baseline for each
%                           neural unit (rows) for each direction (columns).
%       numangles: A number that specifies the number of directions.
%       params: A structure containing the baseline parameters: number of trials,
%               number of units, start and end time, and the direction.
% Output:
%       neural_data: A structure that contains the PSTH, hand position
%                    and hand velocity for each trial for each direction.
    
    for k = 1:numangles
        for i = 1:params.n_units
            for n = 1:params.n_trials
                %create a new structure 'trial1' which returns the same
                %structure as 'trial' with PSTH-baseline, spiketrial -
                %baseline, handPos, handVel.
                neural_data(n,k).PSTH(i,:) = spikedens{i}(k,:) - baseline_spikedens(i,k);
                %neural_data(n,k).spikes(i,:) = trial(n,k).spikes(i,:) - baseline_spikedens(i,k);
                neural_data(n,k).handPos(:,:) = trial(n,k).handPos(:,:);
                neural_data(n,k).baseline = baseline_spikedens;
                for t = 1:length(trial(n,k).spikes(i,:))-1
                     neural_data(n,k).handvel(:,t) = (trial(n,k).handPos(:,t+1) - trial(n,k).handPos(:,t));
                end
                
                % uncomment this in to get acceleration
%                 for t = 1:length(trial(n,k).spikes(i,:))-2
%                      neural_data(n,k).handacc(:,t) = (neural_data(n,k).handvel(:,t+1) - neural_data(n,k).handvel(:,t));
%                 end
            end
        end
    end
end

%% PSTH function

function spikedens = PSTH(trial, params)
% This calculates the Peristimulus time histogram. This is the histograms
% of the times at which neurons fire (which is from 300ms to 500ms: time of movement).
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
%       params: A structure containing the number of trials, number of units,
%               start and end time, and the direction.
% Output:
%       spikedens: For each neural unit, we return the average spike rate
%                  at each ms of movement for each direction.

    for k = params.direction     
        for i = 1:params.n_units
            for n = 1:params.n_trials
                spikecount{i}(n,:) = trial(n,params.direction(k)).spikes(i,params.t_start:params.t_end-1);
            end
            spikedens{i}(k,:) = sum(spikecount{i},1)/params.n_trials;
        end
        
%          Uncomment to plot PSTH
%         figure()
%         for i = 1:params.n_units
%             subplot(params.n_units,1,i)
%             bar(params.t_start:params.t_end-1,spikedens(i,:))
%             hold on
%             axis([params.t_start params.t_end 0 max(max(spikedens))])
%             xlabel('Time (ms)')
%             ylabel(['Neuron = ' num2str(i)])
%         end

        clear spikecount
    end
    
end

%% Baseline function

function baseline_spikedens = baseLine(trial, params)
% This function calculates the baseline (0-300ms of the monkey's non-movement).
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
%       params: A structure containing the number of trials, number of units,
%               start and end time, and the direction.
% Output:
%       baseline_spikedens: Average spike rate of the baseline for each
%                           neural unit (rows) for each direction (columns).
    
    for k = params.direction
        for i = 1:params.n_units
            for n = 1:params.n_trials
                spikecount{i,k}(n,:) = trial(n,params.direction(k)).spikes(i,params.t_start:params.t_end-1);
            end
            spikedens{i,k}(1,:) = sum(spikecount{i,k},1)/params.n_trials;
            baseline_spikedens(i,k) = mean(spikedens{i,k});
        end
        clear spikecount
        clear spikedens
    end
end
