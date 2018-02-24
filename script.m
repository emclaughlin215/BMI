clc
clear all
close all

load monkeydata_training.mat

numangles = 8;
n_units = 5;
angles = [30/180*pi, 70/180*pi, 110/180*pi, 150/180*pi, 190/180*pi,...
230/180*pi, 310/180*pi, 350/180*pi];

% RasterFunct(trial, numangles, n_units)

n_trials = 50;
t_start = 1;
t_end = 300;

PSTH(trial, t_start, t_end, n_units, n_trials, 1:3)

plotXY(trial, angles)