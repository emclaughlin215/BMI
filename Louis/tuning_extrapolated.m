function [preference, std, vector, pref_fit] = tuning_extrapolated(trial, i_id, Plot)
    K=size(trial,2);
    N=size(trial,1);
    I=size(i_id,2);
    preference = zeros(I,K+1);
    std = zeros(I,K+1);
    vector = zeros(I,2);
    vector_norm = zeros(I,2);
        
    for i=i_id
        for n=1:1:N
            for k=1:1:K
                preference(i-i_id(1,1)+1,k)=preference(i-i_id(1,1)+1,k)+sum(trial(n,k).spikes(i,:))/(size(trial(n,k).spikes(i,:),2)*N);
            end
        end
    end
    
    for i=i_id
        for n=1:1:N
            for k=1:1:K
                std(i-i_id(1,1)+1,k)=std(i-i_id(1,1)+1,k)+1/N*(sum(trial(n,k).spikes(i,:))/(size(trial(n,k).spikes(i,:),2))-preference(i-i_id(1,1)+1,k))^2;
            end
        end
    end 
    
    Angles = [30*pi/180:40*pi/180:230*pi/180,310*pi/180,350*pi/180];
     
    ft = fittype('a*sin(x+b)+c','independent','x','dependent','height');
    options = fitoptions(ft);
    options.StartPoint = [1,0,1];
    for i=i_id-i_id(1,1)+1
       pref_fit{i} = fit(Angles',preference(i,1:end-1)',ft,options); 
       preference(i,8:9) = preference(i,7:8);
       preference(i,7) = pref_fit{i}(270*pi/180);
    end
    
    std(:,8:9) = std(:,7:8);
    std(:,7) = 0;
    Angles = 30*pi/180:40*pi/180:350*pi/180;
    Vector = [cos(Angles);sin(Angles)];
    
    for i=i_id-i_id(1,1)+1        
        for k=1:1:K+1
            vector(i,:)=vector(i,:) + preference(i,k)*Vector(:,k)'/sum(preference(i,:));
        end  
        vector_norm(i,:)=vector(i,:)/sqrt(vector(i,1)^2+vector(i,2)^2);
    end
    
    if Plot
        f1=figure(1); set(f1,'name','Preference','numbertitle','off')
        for i=i_id
            subplot(ceil(I/3),3,i-i_id(1,1)+1)
            errorbar(Angles*180/pi,preference(i-i_id(1,1)+1,:),std(i-i_id(1,1)+1,:))
            ylabel(num2str(i))
        end

        f2=figure(2); set(f2,'name','Vector','numbertitle','off')
        for i=i_id
            subplot(ceil(I/3),3,i - i_id(1,1)+1)       
            plot([0,vector_norm(i-i_id(1,1)+1,1)],[0,vector_norm(i-i_id(1,1)+1,2)])
            axis([-1,1,-1,1])
            ylabel(num2str(i-(i-i_id(1,1)+1)))
        end
    end
end