function [decodedPosX, decodedPosY, newParameters] = positionEstimator(trial, Param)
    %Parameters    
    N_iterations = 5;
    speed_std = 0.02 ;
    speed_std2 = 0.3 ;
    t_bin = 20;
    t_planning = 320;
    
    if size(trial.spikes,2)<Param.previous_length
       Param.isfirst = 1;
       Param.decodedPos = trial.startHandPos;
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
        directions_norm = sqrt(Param.direction(:,1).^2+Param.direction(:,2).^2);
        planned_speed = rates'*(Param.direction(:,:)./directions_norm);
        
        %The first step is planification, there is no actual movement
        decodedPosX = trial.startHandPos(1,1);
        decodedPosY = trial.startHandPos(2,1);
        
        newParameters = Param;
        newParameters.Speed_estimate_prev = planned_speed;
        newParameters.particles = planned_speed + randn(Param.N_particles,2)*speed_std2;
        newParameters.decodedPos = [decodedPosX,decodedPosY];
        %We move on to next steps with movement
        newParameters.isfirst = 0;
        newParameters.previous_length = size(trial.spikes,2);
        % Calculate the prefered direction
        angles = [30, 70, 110, 150, 190, 230, 310, 350]/180*pi;
        directions = Param.NET(sum(trial.spikes,2)/size(trial.spikes,2));
        [~, idx] = max(directions);
        newParameters.prefdir = angles(idx);
        newParameters.idx = idx;
    else
        n_mini = 5;
        t_minibin = t_bin/n_mini;
        
        for mini=1:1:n_mini
            %We create a dummy Param structure for the iterations 
            Param_iter = Param;

            %We increment the estimated position
            decodedPosX = Param_iter.decodedPos(1,1) + Param_iter.Speed_estimate_prev(1,1)*t_minibin;
            decodedPosY = Param_iter.decodedPos(1,2) + Param_iter.Speed_estimate_prev(1,2)*t_minibin;

            for iterations=1:1:N_iterations
                %We compute counts (observation)
                counts = sum(trial.spikes(:,end-(n_mini-mini+1)*t_minibin:end-(n_mini-mini)*t_minibin),2);
                %We calculate poisson parameter lambda for each neuron
                Particles_norm = sqrt(sum(Param_iter.particles.^2,2));
                lambda = max(0.0001,Param_iter.baseline(:,1)+Param_iter.direction_sensitivity.*Param_iter.direction*(Param_iter.particles./Particles_norm)'+Param_iter.speed_sensitivity(:,1)*Particles_norm');            
                %Weights calculation (P(observation|state) for each particle)
                weights = zeros(1,Param_iter.N_particles);
                for p=1:1:Param_iter.N_particles
                    weights(1,p) = prod(exp(-lambda(:,p)*t_minibin).*(lambda(:,p)*t_minibin).^counts(:,1)./factorial(counts(:,1)));
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
    %                 plot(Param_iter.particles(:,1),Param_iter.particles(:,2), 'ro')
    %             end
    %             axis([-1 1 -1 1])
    %             pause(0.1)
    %             
                %We add system noise
                Param_iter.particles = randn(Param_iter.N_particles,2)*speed_std + Particles;            
            end
            %After all the iterations, the particle cloud has converged towards
            %the "true" state (i.e. true speed)
            Speed_estimate_prev = mean(Particles);


            %We store parameters for new iteration while adding a -slightly
            %bigger- system noise
            newParameters = Param_iter;
            newParameters.Speed_estimate_prev = correctingSpeed2(Param_iter, Speed_estimate_prev, trial, n_mini, t_minibin, mini);
            newParameters.particles = randn(Param.N_particles,2)*speed_std2 + Particles;
            newParameters.decodedPos = [decodedPosX,decodedPosY];
            newParameters.previous_length = size(trial.spikes,2);
        end
    end
    
end

function newSpeed = correctingSpeed2(Param, v, trial, n_mini, t_minibin, mini)
    longi = [cos(Param.prefdir),sin(Param.prefdir)];
    ortho = [-sin(Param.prefdir),cos(Param.prefdir)];
    x = Param.decodedPos;
    
    k = 0.035;
    k2 = 0.01;
  
%     Magic = [100.118067194997;96.3527537142964;96.3698024072356;95.7070171410154;92.2767943539182;90.8366311233129;97.2950870514154;98.7913981228808];
%     magic = Magic(Param.idx,1);

%       error = abs((x-Param.meanTraj{Param.idx}(1:2,length(trial.spikes)-(n_mini-mini)*t_minibin))*ortho');
%     attractor = magic*longi+trial.startHandPos;
     Magic = [100.118067194997;96.3527537142964;96.3698024072356;95.7070171410154;92.2767943539182;90.8366311233129;97.2950870514154;98.7913981228808];
     magic = Magic(Param.idx,1);
     attractor = magic*longi+trial.startHandPos;
    
    if length(trial.spikes) > length(Param.meanTraj{Param.idx}) 
        error = x*ortho';
    else
%         attractor = Param.meanTraj{Param.idx}(1:2,length(trial.spikes)-(n_mini-mini)*t_minibin);
        sn = sign((x-(Param.meanTraj{Param.idx}(1:2,length(trial.spikes)-(n_mini-mini)*t_minibin)-Param.meanTraj{Param.idx}(1:2,1)+trial.startHandPos))*ortho');
        error = sn.*(x-(Param.meanTraj{Param.idx}(1:2,length(trial.spikes)-(n_mini-mini)*t_minibin)-Param.meanTraj{Param.idx}(1:2,1)+trial.startHandPos))*ortho';
    end
    
    correction = -k2*error*ortho+k*(attractor-x);%/norm(attractor-x);
    
    newSpeed = v+correction;
   
end
    
    
   
