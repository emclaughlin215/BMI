clc
clear all
close all

%load data
load monkeydata_training.mat


numangles = 8; %select number of angles to consider
n_units = 5; %select the number of neurons to consider
angles = [30/180*pi, 70/180*pi, 110/180*pi, 150/180*pi, 190/180*pi,...
230/180*pi, 310/180*pi, 350/180*pi];

% RasterFunct(trial, numangles, n_units)

n_trials = 50; %average over this number of trials
t_start = 300; %start time
t_end = 500; %start time

%calculate the baseline of each neuron
[baseLine] = baseLine(trial);
%calculate the PSTH according to parameters chosen above
PSTH = PSTH(trial, t_start, t_end, 98, 100, 1:8);

for jj = 1:8   
    for j = 1:98
        for i = 1:100
            %create a new structure 'trial1' which returns the same
            %structure as 'trial' with PSTH-baseline, spiketrial -
            %baseline, handPos, handVel.
            trial1(i,jj).PSTH(j,:) = PSTH{j}(jj,:) - baseLine(j,jj);
            trial1(i,jj).spikes(j,:) = trial(i,jj).spikes(j,:) - baseLine(j,jj);
            trial1(i,jj).handPos(:,:) = trial(i,jj).handPos(:,:);
            for ii = 1:length(trial(i,jj).spikes(j,:))-1
                 trial1(i,jj).handvel(:,ii) = (trial(i,jj).handPos(:,ii+1) - trial(i,jj).handPos(:,ii));
            end
% uncomment this in to get acceleration
%             for ii = 1:length(trial(i,jj).spikes(j,:))-2
%                  trial1(i,jj).handacc(:,ii) = (trial1(i,jj).handvel(:,ii+1) - trial1(i,jj).handvel(:,ii));
%             end
        end
    end
end

