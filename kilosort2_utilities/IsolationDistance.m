function [IsolDist,Lratio] = IsolationDistance(Fet, ClusterSpikes, m, toexclude)

% IsolDist = IsolationDistance(Fet, ClusterSpikes, m)
%
% Isolation Distance
% Measure of cluster quality
%
% Inputs:   Fet:           N by D array of feature vectors (N spikes, D dimensional feature space)
%           ClusterSpikes: Index into Fet which lists spikes from the cell whose quality is to be evaluated.
%           m:             squared mahalanobis distances, default is to
%                          calculate them directly
%
% Created by Ken Harris

% find # of spikes in this cluster
ClusterSpikes=ClusterSpikes(:);
toexclude=toexclude(:);
if nargin < 3 || isempty(m)
	nSpikes = size(Fet,1);
else
	nSpikes = size(m,1);
end

nClusterSpikes = length(ClusterSpikes);
InClu = ismember(1:nSpikes, ClusterSpikes);

% check there are enough spikes (but not more than half)
if length(ClusterSpikes) < size(Fet,2) || length(ClusterSpikes)>(nSpikes-length(toexclude))/2
	IsolDist = NaN;
    Lratio = NaN;
	return
end

% mark spikes which are not cluster members
NoiseSpikes = setdiff(1:nSpikes, [ClusterSpikes;toexclude]);

%%%%%%%%%%% compute mahalanobis distances %%%%%%%%%%%%%%%%%%%%%
if nargin < 3 || isempty(m)
	m = mahal(Fet, Fet(ClusterSpikes,:));
end

mCluster = m(ClusterSpikes); % mahal dist of spikes in the cluster
mNoise = m(NoiseSpikes); % mahal dist of all other spikes

% calculate point where mD of other spikes = n of this cell
if nClusterSpikes < (nSpikes-length(toexclude))/2
	[sorted order] = sort(mNoise);
	IsolDist = sorted(nClusterSpikes);
else
	IsolDist = NaN; % If there are more of this cell than every thing else, forget it.
end

df = size(Fet,2);

L = sum(1-chi2cdf(m(NoiseSpikes),df));
Lratio = L/nClusterSpikes;