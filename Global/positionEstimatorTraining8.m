function [Param] = positionEstimatorTraining8(trial_train)
    %We call the filtering function to obtain the data
    trial = filtering_neurons(trial_train, 'None');    
    Param = positionEstimatorTrainingCNN(trial_train);
 
    %Some dimensions for loops
    K = size(trial,2);
    I = size(trial(1).rate,1);
    B = size(trial(1).rate,2);
    
    speed_angle = zeros(K,B);
    
    for k=1:1:K
      for b=1:1:B
        speed_angle(k,b) = atan2(trial(k).speed(2,b),trial(k).speed(1,b));
      end
    end
    
    %Particle filtering parameters
    N_particles = 500;
    
    %Returned values initialization
    direction = zeros(I,2);
    direction_sensitivity = zeros(I,1);
    baseline = trial(1).baseline(:,1);
    
    %We create the point cloud useful for fitting, to obtain rate as a function of speed 
    for i=1:1:I
       Cloud{i}=[0;0];
       for k=1:1:K
            Cloud{i}=[Cloud{i},[speed_angle(k,:);trial(k).rate(i,:)]];
       end
       Cloud{i} = Cloud{i}(:,2:end);
    end    
    
    ft = fittype('exp(a)*cos(x+b)','independent','x','dependent','height');
    options = fitoptions(ft);
    options.StartPoint = [0.1,0];
    for i=1:1:I
        pref_fit{i} = fit(Cloud{i}(1,:)',Cloud{i}(2,:)',ft,options);
        direction(i,:) = [cos(-pref_fit{i}.b),sin(-pref_fit{i}.b)];
        direction_sensitivity(i) = exp(pref_fit{i}.a);
    end
    
    %Retturned parameters
    Param.baseline = baseline;
    Param.direction = direction;
    Param.direction_sensitivity = direction_sensitivity;
    Param.particles = zeros(N_particles,1);
    Param.decodedPos = [0,0];
    Param.isfirst = 1;
    Param.N_particles  = N_particles;
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
    plot(Param.direction(:,1),Param.direction(:,2),'o')
    xlabel('X axis')
    ylabel('Y axis')
    title('Preferred direction') 
    
    neurons_id = 11:13;
    f4 = figure(4);
    f4.Name = 'Neuron fitting';
    angle = -180:1:180;
    for i=neurons_id
        subplot(length(neurons_id),1,i-neurons_id(1)+1)        
        plot(Cloud{i}(1,:)*180/pi,Cloud{i}(2,:),'o',angle,pref_fit{i}(angle*pi/180))
        xlabel('Angle (°)')
        ylabel('Rate (kHz)')
        title(strcat('Neuron ',num2str(i)))    
    end    

end