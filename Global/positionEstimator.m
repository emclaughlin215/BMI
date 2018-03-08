function [decodedPosX, decodedPosY, newParameters] = positionEstimator(trial, Param)
    N_iterations = 10;
    speed_std = 0.1 ;
    speed_std2 = 0.3 ;
    if Param.isfirst
        first_speed_norm = 0.1;       
        
        rates = sum(trial.spikes(:,:),2)/320e-3-Param.baseline;
        directions_norm = sqrt(Param.direction(:,1).^2+Param.direction(:,2).^2);
        planned_speed = first_speed_norm*(sum((Param.direction(:,:)./directions_norm)'*rates));
        
        decodedPosX = 0;
        decodedPosY = 0;
        newParameters = Param;
        newParameters.Particles = planned_speed + randn(Param.N_particles,2)*speed_std;
        newParameters.isfirst = 0;
    else
        Param_iter = Param;
        for iterations=1:1:N_iterations
            counts = sum(trial.spikes(:,end-20:end),2);
            lambda = exp(Param_iter.baseline+Param_iter.direction*Param_iter.particles'+Param_iter.speed_sensitivity.*sqrt(Param_iter.particles(:,1).^2+Param_iter.particles(:,2).^2));
            
            weights = zeros(1,Param_iter.N_particles);
            for p=1:1;Param_iter.N_particles
                weights(1,b) = prod(exp(-lambda(:,p)*20e-3)*(lambda(:,p)*20e-3).^counts(:,1)./factorial(counts(:,1)));
            end
            weights = weights/sum(weights);
            Particles = randsample(Param_iter.particles,Param_iter.N_particles,True,weights);                                    
            Param_iter.particles = randn(Param_iter.N_particles,2)*speed_std + Particles;            
        end
        Speed_estimate = mean(Particles);
        decodedPosX = Param_iter.decodedPos(1,1) + Speed_estimate(1,1)*20e-3;
        decodedPosY = Param_iter.decodedPos(1,2) + Speed_estimate(1,2)*20e-3;
        newParameters = Param_iter;
        newParameters.particles = randn(Param.N_particles,2)*speed_std2 + Particles;
        newParameters.decodedPos = [decodedPosX,decodedPosY];
    end
    
    
   