function [trial_stepping] = data_stepping(trial, N_step)
    N = size(trial, 1);
    K = size(trial, 2);
            
    for n=1:1:N
        for k=1:1:K
            T  =size(trial(n,k).spikes,2);
            t_step=floor(T/N_step);
            for t=1:1:N_step
                trial_stepping(n,k).angles(t) = atan2((trial(n,k).handPos(2,t*t_step)-trial(n,k).handPos(2,(t-1)*t_step+1)),(trial(n,k).handPos(1,t*t_step)-trial(n,k).handPos(1,(t-1)*t_step+1)));
                trial_stepping(n,k).rate(:,t) = sum(trial(n,k).spikes(:,(t-1)*t_step+1:t*t_step),2)/t_step;
            end
        end        
    end

end