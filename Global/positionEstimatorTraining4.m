function [Param] = positionEstimatorTraining4(trial_train)
    %We call the filtering function to obtain the data
    trial = filtering_neurons(trial_train, 'firing_rate');
 
    %Some dimensions for loops
    K = size(trial,2);
    I = size(trial(1).rate,1);
    
    %Particle filtering parameters
    N_particles = 1000;
    
    %Returned values initialization
    direction = zeros(I,2);
    baseline = trial(1).baseline(:,1);
    
    %We create the point cloud useful for fitting, to obtain rate as a function of speed 
    for i=1:1:I
       Cloud{i}=[0;0;0];
       for k=1:1:K
            Cloud{i}=[Cloud{i},[trial(k).speed(1,:);trial(k).speed(2,:);trial(k).rate(i,:)]];
       end
       Cloud{i} = Cloud{i}(:,2:end);
       direction(i,:) = pinv(Cloud{i}(1:2,:)')*Cloud{i}(3,:)';
    end    
    
    %Retturned parameters
    Param.baseline = baseline;
    Param.direction = direction;
    Param.particles = zeros(N_particles,2);
    Param.decodedPos = [0,0];
    Param.isfirst = 1;
    Param.N_particles  =N_particles;
    Param.bool_neurons = trial(1).bool_neurons(:,1);

end