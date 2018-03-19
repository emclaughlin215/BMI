
function [preference, standard_deviation, vector] = tuning(trial, i_id, Plot)
    % It computes (and maybe plots) the Tuning curve for all reaching angles
    % for different neural units for all trials.
    % Input:
    %    trial: A structure that contains the data. 
    %     i_id: A vector that contains all the neural units we want to plot.
    %     Plot: A logical number (1 or 2) that specifies if we want to plot the Tuning curve.
    % Output:
    %           preference: A matrix that is the preferred direction of
    %                       different neural units (rows) for
    %                       different reaching angles (columns).
    %   standard_deviation: A matrix that is the standard deviation across 
    %                       trials for different neural units (rows) for
    %                       different reaching angles (columns).
    %               vector: A matrix that is the x and y direction of the 
    %                       preferred orientation direction for different neural units.
    K = size(trial,2);
    N = size(trial,1);
    I = size(i_id,2);
    preference = zeros(I,K);
    standard_deviation = zeros(I,K);
    vector = zeros(I,2);
    vector_norm = zeros(I,2);
    
    for i = i_id
        for n = 1:1:N
            for k = 1:1:K
                preference(i,k) = preference(i,k)+sum(trial(n,k).rate(i,:))/(size(trial(n,k).rate(i,:),2)*N);
            end
        end
    end
    
    for i = i_id
        for n = 1:1:N
            for k = 1:1:K
                standard_deviation(i,k) = standard_deviation(i,k)+1/N*(sum(trial(n,k).rate(i,:))/(size(trial(n,k).rate(i,:),2))-preference(i,k))^2;
            end
        end
    end
    
    Angles = [30*pi/180:40*pi/180:230*pi/180,310*pi/180,350*pi/180];
    Vector = [cos(Angles);sin(Angles)];
    
    for i = i_id        
        for k = 1:1:K
            vector(i,:) = vector(i,:) + preference(i,k)*Vector(:,k)'/sum(preference(i,:));
        end  
        vector_norm(i,:) = vector(i,:)/sqrt(vector(i,1)^2+vector(i,2)^2);
    end
    
    if Plot
        f1 = figure(1); set(f1,'name','Preference','numbertitle','off')
        for i = i_id
            subplot(ceil(I/3),3,i-i_id(1,1)+1)
            errorbar(Angles*180/pi,preference(i,:),standard_deviation(i,:))
            ylabel(num2str(i))
        end

        f2 = figure(2); set(f2,'name','Vector','numbertitle','off')
        for i = i_id
            subplot(ceil(I/3),3,i-i_id(1,1)+1)       
            plot([0,vector_norm(i,1)],[0,vector_norm(i,2)])
            axis([-1,1,-1,1])
            ylabel(num2str(i))
        end
    end
end
