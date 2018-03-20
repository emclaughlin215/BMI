function [Param] = positionEstimatorTrainingCNN(trial)
    
    rates = zeros(8,800);
    output = zeros(8,800);
    
    
    for k = 1:8
        for j = 1:50
            for i = 1:98 
                rates(i,(k-1)*100+j) = sum(trial(j,k).spikes(i,1:320),2)/320;
                output(k,(k-1)*100+j) = 1;
            end
        end
    end
    
    net = feedforwardnet([10 5 10 5 10 5 10]);
    net = configure(net, rates, output);
    net = init(net);
   [Param.NET, ~] = train(net, rates, output);
    
   %test the accuracy
%    Param.guess = zeros(8,800);
%    for i = 1:800
%     fs = param.NET(rates(:,i));
%     [~, idmax] = max(fs);
%     Param.guess(idmax,i) = 1;
%    end
    
end