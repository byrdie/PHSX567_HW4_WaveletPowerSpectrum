% Octave script to computer wavelet transform for PHSX567

% Load in dataset
data = [csvread('062009.csv');csvread('072009.csv');csvread('082009.csv');csvread('092009.csv')];

% Extract temperature information
temp = data(2:end,7)';

% Generate time information (5 minute intervals)'
time = 0:1:length(temp)-1;

% Compute the mean and variance
mu = mean(temp);
variance = std(temp)^2
