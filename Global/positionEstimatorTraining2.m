function [Param] = positionEstimatorTraining2(trial_train)

    trial = filtering_neurons(trial_train);

    K = size(trial,2);
    I = size(trial(1).rate,1);
    B = size(trial(1).speed,2);
    
    N_particles = 500;
    
    direction = zeros(I,2);
    speed_sensitivity = zeros(I,1);
    baseline = trial(1).baseline(:,1); 
    
    for i=1:1:I
       Cloud{i}=[0;0;0];
       for k=1:1:K
           Cloud{i}=[Cloud{i},[trial(k).speed(1,:);trial(k).speed(2,:);trial(k).rate(i,:)/exp(baseline(i,1))]];
       end
       Cloud{i} = Cloud{i}(:,2:end)
    end
        
    ft = fittype('exp(d1*x+d2*y+s*sqrt(x^2+y^2)','independent',{'x','y'},'dependent','height');
    options = fitoptions(ft);
    options.StartPoint = [1,1,1];
    for i=1:1:I
        pref_fit{i} = fit([Cloud{i}(1,:)',Cloud{i}(2,:)'],Cloud{i}(3,:)',ft,options);
        direction(i,:) = [pref_fit{i}.d1,pref_fit{i}.d2];
        speed_sensitivity(i,1) = pref_fit{i}.s;
    end  
    
    Param.baseline = baseline;
    Param.direction = direction;
    Param.speed_sensitivity = speed_sensitivity;
    Param.particles = zeros(N_particles,2);
    Param.decodedPos = [0;0];
    Param.isfirst = 1;
    Param.N_particles  =N_particles;
    Param.bool_neurons = trial(1).bool_neurons(:,1);

end