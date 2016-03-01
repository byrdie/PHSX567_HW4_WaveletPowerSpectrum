% Octave script to computer wavelet transform for PHSX567

% Load in dataset
data = [csvread('062009.csv');csvread('072009.csv');csvread('082009.csv');csvread('092009.csv')];

% Extract temperature information
temp = data(2:end,7)';

% (2) Calculate wavelet power spectrum and overplot significance contour:
n = length(temp);
dt = 5/(60*24); % sampling period in minutes
t = [0:n-1]*dt;	% time axis
pad = 1;      % pad the time series with zeroes (recommended)
dj = 1/12;    % this will do 12 sub-octaves (semitones in music lingo) per octave.
   % NOTE: I find that oversampling in frequency produces a less "glitchy" reconstruction.
s0 = 2*dt;    % set minimum scale at Nyquist period
j1 = 4.5/dj;    % this says do 7 powers-of-two with dj sub-octaves each
mother = 'Morlet';
Cdelta = 0.776; % Morlet m=6 reconstruction factor.
Psi0_0 = pi^(-1/4); % Morlet mother wavelet amplitude at origin
signif_level = 0.90;

% Some basic signal stats
mu = mean(temp);
variance = std(temp).^2;
lag1 = mean( (temp(2:n)-mu) .* (temp(1:n-1)-mu) ) ./ variance;

printf("The lag-1 autocorrelation is: %f \n", lag1)

% Wavelet transform and significance (Hanford):
[wavlet,period,scale,coi] = wavelet(temp,dt,pad,dj,s0,j1,mother);
frequency = 1./period;
power = (abs(wavlet)).^2 ;        % compute wavelet power spectrum
% Significance levels:
[signif,fft_theor] = wave_signif(variance,dt,scale,0,lag1,signif_level,-1,mother);
sig = (signif').*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig = power ./ sig;         % where ratio > 1, power is significant

figure(2)
hold off
subplot(2,1,1)
imagesc(t,log2(period),log2(power));
%axis(event_times)
hold on
xlabel('time (s)')
ylabel('log_2(period) (s)')
title('Wavelet Power Spectrum')
set(gca,'YDir','reverse')
% significance contour, levels at -99 (fake) and 1 (signif_level)
contour(t,log2(period),sig,[-99, 1],'k');
% cone-of-influence, anything "below" is dubious
plot(t,log2(coi),'k')










