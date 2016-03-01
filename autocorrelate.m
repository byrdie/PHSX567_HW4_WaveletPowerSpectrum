% Function to compute the discrete autocorrection
function R = autocorrelate(obs, mean, variance, k)

	% shift array so we can compare neighboring values
	obs_k = shift(obs, k);
	
	% Eliminate the first value since it will only 
	% complicate things because of the circular shift
	obs = obs(2:end);
	obs_k = obs_k(2:end);
	
	% Save the number of elements for later
	n = length(obs);
	
	% Compute the autocorrelation
	R = (1 ./ (n - k) .* variance) .* sum((obs - mean) .* (obs_k - mean));

endfunction
