% Octave script to computer wavelet transform for PHSX567

% Load in dataset
data = [csvread('062009.csv');csvread('072009.csv');csvread('082009.csv');csvread('092009.csv')];

% Extract temperature information
temp = data(2:end,7)';
clear data;

% (2) Calculate wavelet power spectrum and overplot significance contour:
n = length(temp);
dt = 5/(60*24); % sampling period in minutes
t = [0:n-1]*dt;	% time axis
pad = 1;      % pad the time series with zeroes (recommended)
dj = 1/12;    % this will do 12 sub-octaves (semitones in music lingo) per octave.
   % NOTE: I find that oversampling in frequency produces a less "glitchy" reconstruction.
s0 = 2*dt;    % set minimum scale at Nyquist period
j1 = 13/dj;    % this says do 13 powers-of-two with dj sub-octaves each
mother = 'Morlet';
Cdelta = 0.776; % Morlet m=6 reconstruction factor.
Psi0_0 = pi^(-1/4); % Morlet mother wavelet amplitude at origin
signif_level = 0.90;

% Some basic signal stats
mu = mean(temp);
variance = std(temp).^2;
lag1 = mean( (temp(2:n)-mu) .* (temp(1:n-1)-mu) ) ./ variance;

printf("The lag-1 autocorrelation is: %f \n", lag1)

% Wavelet transform and significance:
[wav,period,scale,coi] = wavelet(temp,dt,pad,dj,s0,j1,mother);
frequency = 1./period;
power = (abs(wav)).^2 ;        % compute wavelet power spectrum
% Significance levels:
[signif,fft_theor] = wave_signif(variance,dt,scale,0,lag1,signif_level,-1,mother);
sig = (signif').*(ones(1,n));  % expand signif --> (J+1)x(N) array
sig = power ./ sig;         % where ratio > 1, power is significant

% Plotting
figure(2)
hold off
imagesc(t,log2(period),log2(power));
%axis(event_times)
hold on
xlabel('time (days)')
ylabel('log_2(period) (days)')
title('Wavelet Power Spectrum')
set(gca,'YDir','reverse')
set (gca, "dataaspectratio", [16 3 1])
% significance contour, levels at -99 (fake) and 1 (signif_level)
contour(t,log2(period),sig,[-99, 1],'k',"linewidth", 2);
% cone-of-influence, anything "below" is dubious
plot(t,log2(coi),'k',"linewidth", 2)

% Horizontal lines marking harmonics
a = log2(1)*ones(1,n);
b = log2(1/2)*ones(1,n);
c = log2(1/3)*ones(1,n);
d = log2(1/4)*ones(1,n);
plot(t, a, 'k',"linewidth", 3);
plot(t, b, 'k',"linewidth", 3);
plot(t, c, 'k',"linewidth", 3);
plot(t, d, 'k',"linewidth", 3);

% Filter and reconstruct the time series
lows = find(sig < 0.4);
w_filt = wav;
w_filt(lows) = 0;
filt = (dj*sqrt(dt)/Cdelta/Psi0_0) .* sum( real(w_filt) ./ (sqrt(scale')*ones(1,n)), 1 ) ;

% plot the reconstructed time series
figure(3);
plot(t, filt, 'b', 'linewidth', 2.0)
title('Reconstructed time series')
xlabel('time (days)')
ylabel('Temperature')






