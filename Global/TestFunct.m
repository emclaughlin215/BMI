clear all
close all
clc

load monkeydata_training

numNeurons =98;
numParticles = 500;
numbins = 5;
duration = 100; %ms
cutoffweight = .5;
sigmaV = 0.01;

neuronProps = zeros(numNeurons,4);
neuronProps(:,1) = randi(200,numNeurons,1); %baseline
neuronProps(:,2) = randi(200,numNeurons,1); %sensitivity
neuronProps(:,3:4) = [randi(360,numNeurons,1), randi(360,numNeurons,1)]; %preferred angle

% neuronProps(:,1) = param.baseline; %baseline
% neuronProps(:,2) = param.sensitivity; %sensitivity
% neuronProps(:,3:4) = param.direction; %preferred angle

V(:,:) = [2*rand(1,numParticles)-1; 2*rand(1,numParticles)-1]';
figure(1)
scatter(V(:,1),V(:,2),'r')
xlabel('x velocity')
ylabel('y velocity')
grid on 
pause

for n = 1:100/20

    V(:,:) = [2*rand(1,numParticles)-1; 2*rand(1,numParticles)-1]';
    
    for m = 1:10

        lambda(:,:) = max(neuronProps(:,1) + neuronProps(:,2).*(neuronProps(:,3:4)*V'),0);

        for i = 1:numNeurons
            y = sum(trial(50,1).spikes(i,20*(n-1)+1:20*n));
            poissondist(i,:) = exp(-lambda(i,:)*0.02).*lambda(i,:)*0.02 ^y / factorial(y);
        end

        weights = sum(poissondist,1);

        NewParts = [];
        for j = 1:numParticles
            if (weights(j)*numParticles) < cutoffweight
                NewParts = [NewParts;j];
                weights(j) = 0;
            end
        end

        %assign the new particles to ones with higher weights
        %by copying values 
        for j = 1:numParticles
            %for each new particle, copy a particle with a high weight
            % weighted random sample of indexes by their normalized weight:
            particle_to_copy = randsample(1:numParticles, 1, true, weights);
            VNEW(j,:) = V(particle_to_copy,:) + [sigmaV*randn(1,1), sigmaV*randn(1,1)];
        end

        V = VNEW;

        figure(1)
        scatter(V(:,1),V(:,2),'r')
        xlabel('x velocity')
        ylabel('y velocity')
        grid on 
        
        pause(0.5)

    end
    
    figure(2)
    scatter(V(:,1),V(:,2),'b')
    xlabel('x velocity')
    ylabel('y velocity')
    grid on
    
    pause
    
end