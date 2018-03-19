function [Param] = positionEstimatorTraining(trial_train)
    %We call the filtering function to obtain the data
    trial = filtering_neurons(trial_train, 'None');
 
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
    N_particles = 800;
    
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