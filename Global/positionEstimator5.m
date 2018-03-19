function [decodedPosX, decodedPosY, newParameters] = positionEstimator(trial, Param)
    %Parameters    
    N_iterations = 15;
    speed_std = 0.02 ;
    speed_std2 = 0.3 ;
    t_bin = 20;
    t_planning = 320;
    
    if size(trial.spikes,2)<Param.previous_length
       Param.isifirst = 1;
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
        first_speed_norm = 0.1;       
        
        %We obtain the rates and normalized directions, and make a weighted
        %sum of the latter 
        rates = sum(trial.spikes(:,:),2)/t_planning-Param.baseline;
        directions_norm = sqrt(Param.direction(:,1).^2+Param.direction(:,2).^2);
        planned_speed = first_speed_norm*rates'*(Param.direction(:,:)./directions_norm);
        
        %The first step is planification, there is no actual movement
        decodedPosX = trial.startHandPos(1,1);
        decodedPosY = trial.startHandPos(2,1);
        
        newParameters = Param;
        newParameters.Speed_estimate_prev = planned_speed;
        newParameters.particles = planned_speed + randn(Param.N_particles,2)*speed_std2;
        %We move on to next steps with movement
        newParameters.isfirst = 0;
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
            Particles_norm = sqrt(sum(Param_iter.particles.^2,2));
            lambda = max(0,Param_iter.baseline(:,1)+Param_iter.direction_sensitivity.*Param_iter.direction*(Param_iter.particles./Particles_norm)'+Param_iter.speed_sensitivity(:,1)*Particles_norm');            
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
            f3 = figure(3);
            f3.Name = 'Speed particles population';
            if iterations == N_iterations
                plot(Param_iter.particles(:,1),Param_iter.particles(:,2), 'ro')
            end
            axis([-1 1 -1 1])
            pause(0.1)
            
            %We add system noise
            Param_iter.particles = randn(Param_iter.N_particles,2)*speed_std + Particles;            
        end
        %After all the iterations, the particle cloud has converged towards
        %the "true" state (i.e. true speed)
        Speed_estimate_prev = mean(Particles);
        
        %We store parameters for new iteration while adding a -slightly
        %bigger- system noise
        newParameters = Param_iter;
        newParameters.Speed_estimate_prev = Speed_estimate_prev;
        newParameters.particles = randn(Param.N_particles,2)*speed_std2 + Particles;
        newParameters.decodedPos = [decodedPosX,decodedPosY];
        newParameters.previous_length = size(trial.spikes,2);
    end
    
end
    
    
   