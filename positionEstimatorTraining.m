function [Param] = positionEstimatorTraining(trial_train)
    %We call the filtering function to obtain the data
    trial = filtering_neurons(trial_train, 'None');
    Param = positionEstimatorTrainingCNN(trial_train);
 
    %Some dimensions for loops
    K = size(trial,2);
    I = size(trial(1).rate,1);
    B = size(trial(1).rate,2);
    
    speed_angle = zeros(K,B);
    speed_norm = zeros(K,B);
    
    for k=1:1:K
      for b=1:1:B
        speed_angle(k,b) = atan2(trial(k).speed(2,b),trial(k).speed(1,b));
        speed_norm(k,b) = sqrt(trial(k).speed(2,b)^2+trial(k).speed(1,b)^2);
      end
    end
    
    %Particle filtering parameters
    N_particles = 100;
    
    %Returned values initialization
    direction = zeros(I,2);
    speed_sensitivity = zeros(I,1);
    direction_sensitivity = zeros(I,1);
    baseline = trial(1).baseline(:,1);
    
    %We create the point cloud useful for fitting, to obtain rate as a function of speed 
    for i=1:1:I
       Cloud{i}=[0;0;0];
       for k=1:1:K
            Cloud{i}=[Cloud{i},[speed_angle(k,:);speed_norm(k,:);trial(k).rate(i,:)]];
       end
       Cloud{i} = Cloud{i}(:,2:end);
    end    
    
    ft = fittype('exp(a)*cos(x+b)','independent','x','dependent','height');
    options = fitoptions(ft);
    options.StartPoint = [0.1,0];
    for i=1:1:I
        pref_fit{i} = fit(Cloud{i}(1,:)',Cloud{i}(3,:)',ft,options);
        direction(i,:) = [cos(-pref_fit{i}.b),sin(-pref_fit{i}.b)];
        direction_sensitivity(i) = exp(pref_fit{i}.a);
        %speed_sensitivity(i,1) = pinv(Cloud{i}(2,:)')*Cloud{i}(3,:)';
        speed_sensitivity(i,1) = pinv(Cloud{i}(2,:)')*(Cloud{i}(3,:)'-pref_fit{i}(Cloud{i}(1,:)'));
    end
    
    %Retturned parameters
    Param.baseline = baseline;
    Param.direction = direction;
    Param.direction_sensitivity = direction_sensitivity;
    Param.speed_sensitivity = speed_sensitivity;
    Param.particles = zeros(N_particles,2);
    Param.decodedPos = [0,0];
    Param.isfirst = 1;
    Param.N_particles  =N_particles;
    Param.bool_neurons = trial(1).bool_neurons(:,1);
    Param.previous_length = 0;
    
    %Plot
    f2 = figure(2);
    f2.Name = 'Neurons characteristics';
    subplot(2,2,1)
    histogram(Param.baseline)
    ylabel('Neuron count')
    xlabel('Baseline (kHz)')
    title('Baselines')
    subplot(2,2,2)
    histogram(Param.direction_sensitivity)
    ylabel('Neuron count')
    xlabel('Direction sensitivity (kHz)')
    title('Direction sensitivity')
    subplot(2,2,3)
    histogram(Param.speed_sensitivity)
    ylabel('Neuron count')
    xlabel('Speed sensitivity (kHz.ms/m)')
    title('Speed sensitivity')
    subplot(2,2,4)
    plot(Param.direction(:,1),Param.direction(:,2),'o')
    xlabel('X axis')
    ylabel('Y axis')
    title('Preferred direction') 
    
    neurons_id = 11:13;
    f4 = figure(4);
    f4.Name = 'Neuron fitting';
    angle = -180:1:180;
    speed = 0:0.01:0.9;
    for i=neurons_id
        subplot(length(neurons_id),2,2*(i-neurons_id(1))+1)        
        plot(Cloud{i}(1,:)*180/pi,Cloud{i}(3,:),'o',angle,pref_fit{i}(angle*pi/180))
        xlabel('Angle (°)')
        ylabel('Rate (kHz)')
        title(strcat('Neuron ',num2str(i)))
        subplot(length(neurons_id),2,2*(i-neurons_id(1))+2)
        plot(Cloud{i}(2,:),Cloud{i}(3,:),'o',speed,speed_sensitivity(i,1)*speed)
        xlabel('Speed norm (m/ms)')
        ylabel('Rate (kHz)')
        title(strcat('Neuron ',num2str(i)))       
    end    

end

function [Param] = positionEstimatorTrainingCNN(trial)
    
    rates = zeros(8,800);
    output = zeros(8,800);
    Param.meanTraj = {};
    J = size(trial,1);
    K = size(trial,2);
    I = size(trial(1,1).spikes,1);
    
    for k = 1:K
        Param.meanTraj{k} = zeros(2,550);
        for j = 1:J
            for i = 1:I 
                rates(i,(k-1)*100+j) = sum(trial(j,k).spikes(i,1:320),2)/320;
                output(k,(k-1)*100+j) = 1;
            end
            deltaMeanTraj = trial(j,k).handPos(1:2,1:550)/50;
            Param.meanTraj{k} = Param.meanTraj{k} + deltaMeanTraj;
        end
    end
    
    net = feedforwardnet([10 5 10 5 10 5 10]);
    net = configure(net, rates, output);
    net = init(net);
   [Param.NET, ~] = train(net, rates, output);
   
    
   %test the accuracy
%    Param.guess = zeros(8,800);
%    for i = 1:800
%     fs = param.NET(rates(:,i));
%     [~, idmax] = max(fs);
%     Param.guess(idmax,i) = 1;
%    end
    
end


function filtered_trial = filtering_neurons(trial, type)
% It filters the neurons that we want to use for our analysis.
% We decide if we take into account the neuron or not based on if their firing
% rate is higher compared to the baseline firing rate (before movement).
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles
%       type: A string that can be the type of filtering we want to do (FF, Partial)
% Output:
%       filtered_trial: A structure that contains baseline, rate and speed (divided into bins) for each orientation. 

    bins = 20; % number of divisions we want
    neural_data = getNeuronData(trial, bins);
    baseline = neural_data{1}.baseline;
    
    num_angles = size(trial,2); % number angles
    num_trials = size(trial,1); % number trials

    for i = 1:size(baseline,1) % neuron
        rate_movement_orients = [];
        for k = 1:num_angles % orientation
            % Average firing rate for the movement phase for this neuron across all trials
            rate_movement_orients(k,:) = neural_data{k}.PSTH(i,:);
        end

        % Variance of firing rate for this neuron across all orientations
        var_neuron_orient(i,:) = var(rate_movement_orients);
        
        % Average variance across all orientations and across all trials for this neural unit
        var_movement_neuron(i,1) = mean(var_neuron_orient(i,:));
        
        % Average firing rate across all directions and across all trials for this neural unit
        avg_baseline_neuron(i,1) = mean(baseline(i,:)); % 0 - 300ms
        avg_movement_neuron(i,1) = mean(mean(rate_movement_orients)); % 300 - 500ms

        % Different measures
        FF(i,1) = var_movement_neuron(i)/avg_baseline_neuron(i);
        firing_rate_ratio(i,1) = avg_movement_neuron(i)/avg_baseline_neuron(i);
    end
    
    % Decide if we want to get rid of the neuron
    switch type
        case 'FF'
            outliers_idx = isoutlier(FF);
            for neuron = 1:length(outliers_idx)
                if outliers_idx(neuron) == 1 && outliers_idx(neuron) < mean(FF)
                    bool_neurons(neuron,1) = 0; % we remove
                else
                    bool_neurons(neuron,1) = 1; % we keep
                end
            end
        case 'firing_rate'
            for neuron = 1:size(baseline,1)
                if firing_rate_ratio(neuron) < 0 % we get rid of neuron
                    bool_neurons(neuron,1) = 0; % we remove
                else
                    bool_neurons(neuron,1) = 1; % we keep
                end
            end
        case 'None'
            for neuron = 1:size(baseline,1)
               bool_neurons(neuron,1) = 1; 
            end
    end
    bool_neurons = logical(bool_neurons);
    
    % Remove rate and baseline of this neuron for all
    % orientations and all trials
    for k = 1:num_angles
        neural_data{k}.PSTH(~bool_neurons,:) = [];
    end
    avg_baseline_neuron(~bool_neurons,:) = [];    
    
    % Average PSTH across all trials for each direction for each neuron
    for k = 1:num_angles
        rate{k} = neural_data{k}.PSTH;
        avg_velocity{k} = neural_data{k}.handVel;
    end
    
    [rate_split, velocity_split] = data_stepping(trial, bins, rate, avg_velocity);
 
    % Output
    for k = 1:size(baseline,2)
        filtered_trial(k).baseline = avg_baseline_neuron;
        filtered_trial(k).rate = rate_split{k};
        filtered_trial(k).speed = velocity_split{k};
        filtered_trial(k).bool_neurons = bool_neurons;
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

function neural_data = getNeuronData(trial,bins)
% Calculates hand velocity, hand position, spikes without baseline, 
% PSTH and baseline for every trial in each direction.
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
%       bins: Number of divisions of data.
% Output:
%       neural_data: A structure that contains the PSTH, hand position
%                    and hand velocity for each trial for each direction.

    numangles = 8; % select number of angles to consider
    
    %% Baseline
    % Parameters for baseline function
    params_baseline.n_trials = size(trial,1); % number trials
    params_baseline.n_units = size(trial(1,1).spikes,1); % number units
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
    neural_data = handVelocity(trial, spikedens, baseline_spikedens, numangles, params_baseline, params_PSTH, bins);
end

function neural_data = handVelocity(trial, spikedens, baseline_spikedens, numangles, params, params_PSTH, bins)
% This calculates the hand velocity and hand position averaged across all
% trials for each orientation.
% Input:
%       trial: A structure that contains the data (100 trials across 8 angles)
%       spikedens: For each neural unit, we return the average spike rate
%                  at each ms of movement for each direction.
%       baseline_spikedens: Average spike rate of the baseline for each
%                           neural unit (rows) for each direction (columns).
%       numangles: A number that specifies the number of directions.
%       params: A structure containing the baseline parameters: number of trials,
%               number of units, start and end time, and the direction.
%       bins: Number of divisions of data.
% Output:
%       neural_data: A structure that contains the PSTH, baseline, hand position
%                    and hand velocity across all trials for each direction.
    
    % Find the longest handPos for each orientation
    time_movement = params_PSTH.t_start:1:params_PSTH.t_end;
    
    % Implement delay by taking the next bin_count ms of velocity
    T = time_movement(end)-time_movement(1); % time
    bin_count = floor(T/bins);
    time_begin = time_movement(1) + bin_count;
    time_end = time_movement(end) + bin_count;
    time_movement = time_begin:time_end;
    
    for k = 1:numangles
        max_length_pos(k) = -inf; % initialize
        for n = 1:params.n_trials
            if length(trial(n,k).handPos(1,time_movement)) > max_length_pos(k)
                max_length_pos(k) = length(trial(n,k).handPos(1,time_movement));
            end
        end
    end
    
    % Average handPos and handVel across all trials
    for k = 1:numangles
        % Position
        handPos_x = NaN(params.n_trials,max_length_pos(k));
        handPos_y = NaN(params.n_trials,max_length_pos(k));
        handPos_z = NaN(params.n_trials,max_length_pos(k));
        % Velocity
        handVel_x = NaN(params.n_trials,max_length_pos(k)-1);
        handVel_y = NaN(params.n_trials,max_length_pos(k)-1);
        handVel_z = NaN(params.n_trials,max_length_pos(k)-1);
        for n = 1:params.n_trials
            current_length = length(trial(n,k).handPos(1,time_movement));
            handVel = diff(trial(n,k).handPos(:,time_movement),1,2); % obtain hand velocity
            
            handPos_x(n,1:current_length) = trial(n,k).handPos(1,time_movement);
            handPos_y(n,1:current_length) = trial(n,k).handPos(2,time_movement);
            handPos_z(n,1:current_length) = trial(n,k).handPos(3,time_movement);
            
            handVel_x(n,1:current_length-1) = handVel(1,:);
            handVel_y(n,1:current_length-1) = handVel(2,:);
            handVel_z(n,1:current_length-1) = handVel(3,:);
        end
        neural_data{k}.handPos = [nanmean(handPos_x); nanmean(handPos_y); nanmean(handPos_z)];
        neural_data{k}.handVel = [nanmean(handVel_x); nanmean(handVel_y); nanmean(handVel_z)];
    end
    
    for k = 1:numangles
        for i = 1:params.n_units
            neural_data{k}.PSTH(i,:) = spikedens{i}(k,:) - baseline_spikedens(i,k);
            % Remove negative values and set to 0
            %idx = find(neural_data{k}.PSTH(i,:) < 0);
            %neural_data{k}.PSTH(i,idx) = 0;
        end
        neural_data{k}.baseline = baseline_spikedens;
    end
end

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
