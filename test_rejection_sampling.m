function [distr,overlap,srate] = test_rejection_sampling(samples,target,options)
% [distr,overlap,srate] = TEST_REJECTION_SAMPLING(samples,target,options)
%
% Reproduces a target distribution by von Neumann rejection sampling from a
% distribution of samples and computes the overlap with the target
%
% Input:
%
% samples   assumed distribution of original samples, vector (N,1) of data
%           points, g(r)
% target    target distribution, vector (N,1) of data points, defaults to
%           samples, f(r)
% options   computation options, struct with fields
%           .trials     number of Monte Carlo trials
%
% Output:
%
% distr     output distribution, vector (N,1)
% overlap   overlap of output and target distribution
% srate     success rate, number of successful trials divided by total
%           number of trials

% Initialize empty output
distr = [];
overlap = [];

% Normalize samples distribution and compute cumulative sum
samples = samples/sum(samples);
cum_samples = cumsum(samples);

% initialize unset input
if ~exist('target','var') || isempty(target)
    target = samples;
end

if ~exist('options','var') || ~isfield(options,'trials') || isempty(options.trials)
    options.trials = 1e6;
end

% Check input consistency, return if inconsistent
[N,~] = size(samples);
[N2,~] = size(target);
if N ~= N2 % error handling
    fprintf(2,'Inconsistent length of samples and target vectors (%i,%i). Aborting\n',N,N2);
    return
end

% initialize zero output distribution
distr= zeros(size(target));

% Monte Carlo loop
success = 0;
for trial = 1:options.trials
    [~,point] = min(abs(cum_samples - rand));
    ACCEPT = true;
    if ACCEPT
        distr(point) = distr(point) + 1;
        success = success + 1;
    end
end

% normalize output distribution
distr = distr/sum(distr);

% compute success rate
srate = success/options.trials;

% compute overlap
overlap = sum(min([distr;target]));