function [decodedPosX, decodedPosY, newParameters] = positionEstimator8(trial, Param)
    %Parameters    
    N_iterations = 10;
    angle_std = 1*pi/180 ;
    angle_std2 = 10*pi/180 ;
    t_bin = 20;
    t_planning = 320;
    t_end = 1000;
    t_vmu = 410;
    speed_std = 100;
    speed_ratio = 103.9008/0.9993;
    
    time = t_planning:t_bin:t_end;
    speed = normpdf(time,t_vmu,speed_std)*speed_ratio; 
    
    if size(trial.spikes,2)<Param.previous_length
       Param.isfirst = 1;
       Param.decodedPos = [0,0];
    end
    
    %Neuron filtering
    N = size(trial,1);
    K = size(trial,2);
    for n=1:1:N
       for k=1:1:K
          trial(n,k).spikes = trial(n,k).spikes(Param.bool_neurons,:); 
       end
    end
    
    if Param.isfirst
        %For first estimate we use poplation vector as the expected value
        %for a Gaussian repartition of particles
                
        %We obtain the rates and normalized directions, and make a weighted
        %sum of the latter 
        rates = sum(trial.spikes(:,:),2)/t_planning-Param.baseline;
        planned_speed = speed(1,1)*rates'*Param.direction(:,:);
        
        %The first step is planification, there is no actual movement
        decodedPosX = trial.startHandPos(1,1);
        decodedPosY = trial.startHandPos(2,1);
        
        newParameters = Param;
        newParameters.Speed_estimate_prev = planned_speed;
        newParameters.particles = atan2(planned_speed(1,2),planned_speed(1,1)) + randn(Param.N_particles,1)*angle_std2;
        %We move on to next steps with movement
        newParameters.isfirst = 0;
        newParameters.previous_length = size(trial.spikes,2);
    else
        %We create a dummy Param structure for the iterations 
        Param_iter = Param;
        
        %We increment the estimated position
        decodedPosX = Param_iter.decodedPos(1,1) + Param_iter.Speed_estimate_prev(1,1)*t_bin;
        decodedPosY = Param_iter.decodedPos(1,2) + Param_iter.Speed_estimate_prev(1,2)*t_bin;
        
        for iterations=1:1:N_iterations
            %We compute counts (observation)
            counts = sum(trial.spikes(:,end-t_bin:end),2);
            %We calculate poisson parameter lambda for each neuron
            lambda = max(0.0001,Param_iter.baseline(:,1)+Param_iter.direction_sensitivity.*Param_iter.direction*[cos(Param_iter.particles),sin(Param_iter.particles)]');            
            %Weights calculation (P(observation|state) for each particle)
            weights = zeros(1,Param_iter.N_particles);
            for p=1:1:Param_iter.N_particles
                weights(1,p) = prod(exp(-lambda(:,p)*t_bin).*(lambda(:,p)*t_bin).^counts(:,1)./factorial(counts(:,1)));
            end
            %We resample particles according to the weights: "survival of
            %the fittest"
            weights = weights/sum(weights);
            PartIdx = randsample(1:length(Param_iter.particles),Param_iter.N_particles,true,weights);
            Particles = Param_iter.particles(PartIdx,:);
            
            %This plot helps to show whats happening.
%             f3 = figure(3);
%             f3.Name = 'Speed particles population';
%             if iterations == N_iterations
%                 plot(cos(Param_iter.particles),sin(Param_iter.particles), 'ro')
%             end
%             axis([-1 1 -1 1])
%             pause(0.1)
            
            %We add system noise
            Param_iter.particles = randn(Param_iter.N_particles,1)*angle_std + Particles;            
        end
        %After all the iterations, the particle cloud has converged towards
        %the "true" state (i.e. true speed)
        Speed_estimate_prev = speed(1,(size(trial.spikes,2)-300)/t_bin)*[cos(mean(Particles)),sin(mean(Particles))];
        
        %We store parameters for new iteration while adding a -slightly
        %bigger- system noise
        newParameters = Param_iter;
        newParameters.Speed_estimate_prev = Speed_estimate_prev;
        newParameters.particles = randn(Param.N_particles,1)*angle_std2 + Particles;
        newParameters.decodedPos = [decodedPosX,decodedPosY];
        newParameters.previous_length = size(trial.spikes,2);
    end
    
end
    
    
   