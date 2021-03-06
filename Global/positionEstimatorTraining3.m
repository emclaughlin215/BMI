function [Param] = positionEstimatorTraining3(trial_train)
    %We call the filtering function to obtain the data
    trial = filtering_neurons(trial_train, 'firing_rate');
 
    %Some dimensions for loops
    K = size(trial,2);
    I = size(trial(1).rate,1);
    
    %Particle filtering parameters
    N_particles = 300;
    
    %Returned values initialization
    direction = zeros(I,2);
    speed_sensitivity = zeros(I,1);
    baseline = trial(1).baseline(:,1);
    baseline_other = zeros(I,1);
    
    %We create the point cloud useful for fitting, to obtain rate as a function of speed 
    for i=1:1:I
       Cloud{i}=[0;0;0];
       for k=1:1:K
            Cloud{i}=[Cloud{i},[trial(k).speed(1,:);trial(k).speed(2,:);(trial(k).rate(i,:)+baseline(i,1))]];
       end
       Cloud{i} = Cloud{i}(:,2:end);
    end
    
    %Function fitting, see documentation. 
    ft = fittype('exp(b+d1*x+d2*y+s*sqrt(x^2+y^2))','independent',{'x','y'},'dependent','height');
    options = fitoptions(ft);
    options.StartPoint = [0.5,1,1,1];
    for i=1:1:I
        pref_fit{i} = fit([Cloud{i}(1,:)',Cloud{i}(2,:)'],Cloud{i}(3,:)',ft,options);
        direction(i,:) = [pref_fit{i}.d1,pref_fit{i}.d2];
        speed_sensitivity(i,1) = pref_fit{i}.s;
        baseline_other(i,1) = pref_fit{i}.b;
    end  
    
    %Retturned parameters
    Param.baseline = baseline_other;
    Param.direction = direction;
    Param.speed_sensitivity = speed_sensitivity;
    Param.particles = zeros(N_particles,2);
    Param.decodedPos = [0,0];
    Param.isfirst = 1;
    Param.N_particles  =N_particles;
    Param.bool_neurons = trial(1).bool_neurons(:,1);
    
     %Plot
    f2 = figure(2);
    f2.Name = 'Neurons characteristics';
    subplot(2,2,1)
    histogram(Param.baseline)
    ylabel('Neuron count')
    xlabel('Baseline (kHz)')
    title('Baselines')
    subplot(2,2,2)
    histogram(Param.speed_sensitivity)
    ylabel('Neuron count')
    xlabel('Speed sensitivity (kHz.ms/m)')
    title('Speed sensitivity')
    subplot(2,2,3)
    plot(Param.direction(:,1),Param.direction(:,2),'o')
    xlabel('X axis')
    ylabel('Y axis')
    title('Preferred direction')
    
    neurons_id = 11:13;
    f4 = figure(4);
    f4.Name = 'Neuron fitting';
    vx=-0.9:0.01:0.9;
    vy=-0.9:0.01:0.9;
    [meshVx,meshVy] = meshgrid(vx,vy);
    for i=neurons_id
        subplot(length(neurons_id),1,i-neurons_id(1)+1)        
        plot3(Cloud{i}(1,:),Cloud{i}(2,:),Cloud{i}(3,:),'o')
        hold on
        meshFit = pref_fit{i}(meshVx,meshVy);
        surf(meshVx,meshVy,meshFit);
        xlabel('X speed (m/ms)')
        ylabel('Y speed (m/ms)')
        zlabel('Rates (kHz)')
        title(strcat('Neuron ',num2str(i)))   
    end    

end