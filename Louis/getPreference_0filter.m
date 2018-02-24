function [Cloud, pref_fit, vector, vector_norm] = getPreference_0filter(trial,N_step,cloud_filter,Plot)
    [trial_stepping] = data_stepping(trial, N_step);
    N = size(trial,1);
    K = size(trial,2);
    I = size(trial(1,1).spikes,1);
    vector = zeros(I,2);
    vector_norm = zeros(I,2);
    
    for i=1:1:I
        Cloud{i} = [0;0];
        for n=1;1:N
            for k=1:1:K
                if trial_stepping(n,k).rate(i,:)~=0
                    Cloud{i}=[Cloud{i},[trial_stepping(n,k).angles(1,:);trial_stepping(n,k).rate(i,:)]];
                end
            end
        end
        Cloud{i} = Cloud{i}(:,2:end);
    end
    
    ft = fittype('a*sin(x+b)+c','independent','x','dependent','height');
    options = fitoptions(ft);
    options.StartPoint = [1,0,1];
    
    for i=1:1:I
        if size(Cloud{i},2)>max(cloud_filter,2)
            pref_fit{i} = fit(Cloud{i}(1,:)',Cloud{i}(2,:)',ft,options);
        else
            pref_fit{i} = @(x)0;
        end
    end
    
    for i=1:1:I
        f_cos = @(x)pref_fit{i}(x)*cos(x);
        f_sin = @(x)pref_fit{i}(x)*sin(x);
        vector(i,:) = [integral(f_cos,-pi,pi,'ArrayValued',true),integral(f_sin,-pi,pi,'ArrayValued',true)];
        if norm(vector(i,:))~=0
            vector_norm(i,:) = vector(i,:)/norm(vector(i,:));
        else
            vector_norm(i,:) = vector(i,:);
        end
    end
        
    if Plot
        Angles = -pi:1*pi/360:pi;
        f1=figure(1); set(f1,'name','Preference function','numbertitle','off');
        for i=1:1:I
            subplot(ceil(I/3),3,i)
            plot(Angles*180/pi,pref_fit{i}(Angles),Cloud{i}(1,:)*180/pi,Cloud{i}(2,:),'o')
            axis([-180,180,-Inf,Inf])
            ylabel(num2str(i))
            if size(Cloud{i},2)<=cloud_filter
                fill([-200,-200,200,200],[-1,1,1,-1],'b')
            end
        end
        f2=figure(2); set(f2,'name','Vector','numbertitle','off');
        for i=1:1:I
            subplot(ceil(I/3),3,i)
            plot([0,vector_norm(i,1)],[0,vector_norm(i,2)])
            axis([-1,1,-1,1])
            ylabel(num2str(i))
            if size(Cloud{i},2)<=cloud_filter
                fill([-1,-1,1,1],[-1,1,1,-1],'b')
            end
        end
    end    

end