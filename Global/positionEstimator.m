function [decodedPosX, decodedPosY, newParameters] = positionEstimator(trial, Param)
    %Parameters    
    N_iterations = 1;
    speed_std = 0.01 ;
    speed_std2 = 0.05 ;
    t_bin = 20;
    t_planning = 320;
    
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
        decodedPosX = 0;
        decodedPosY = 0;
        newParameters = Param;
        newParameters.particles = planned_speed + randn(Param.N_particles,2)*speed_std2;
        %We move on to next steps with movement
        newParameters.isfirst = 0;
    else
        %We create a dummy Param structure for the iterations 
        Param_iter = Param;
        for iterations=1:1:N_iterations
            %We compute counts (observation)
            counts = sum(trial.spikes(:,end-t_bin:end),2);
            %We calculate poisson parameter lambda for each neuron
            lambda = exp(Param_iter.baseline+Param_iter.direction*Param_iter.particles'+Param_iter.speed_sensitivity*sqrt(Param_iter.particles(:,1).^2+Param_iter.particles(:,2).^2)');
            
            %Weights calculation (P(observation|state) for each particle)
            weights = zeros(1,Param_iter.N_particles);
            for p=1:1:Param_iter.N_particles
                weights(1,p) = prod(exp(-lambda(:,p)*t_bin).*(lambda(:,p)*t_bin).^counts(:,1)./factorial(counts(:,1)));
            end
            %We resample particles according to the weights: "survival of
            %the fittest"
            weights = weights/sum(weights);
            Particles = datasample(Param_iter.particles,Param_iter.N_particles,1,'Replace',true,'Weights',weights'); 
            %We add system noise
            Param_iter.particles = randn(Param_iter.N_particles,2)*speed_std + Particles;            
        end
        %After all the iterations, the particle cloud has converged towards
        %the "true" state (i.e. true speed)
        Speed_estimate = mean(Particles);
        %We increment the estimated position
        decodedPosX = Param_iter.decodedPos(1,1) + Speed_estimate(1,1)*t_bin;
        decodedPosY = Param_iter.decodedPos(1,2) + Speed_estimate(1,2)*t_bin;
        %We store parameters for new iteration while adding a -slightly
        %bigger- system noise
        newParameters = Param_iter;
        newParameters.particles = randn(Param.N_particles,2)*speed_std2 + Particles;
        newParameters.decodedPos = [decodedPosX,decodedPosY];
    end
    
end
    
    
   