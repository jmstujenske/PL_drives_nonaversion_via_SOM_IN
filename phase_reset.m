function [phase_reset_strength,phases,raw]=phase_reset(signal,signal_times,event_times,before,after,freq_range,fs,order)
%[phase_reset_strength,phases,out]=phase_reset(signal,signal_times,event_times,before,after,freq_range,fs)
%
%Requires circular statistics toolbox!
%
%Calculates theta-phase reset strength
%IN:
%signal - raw signal; can be any number of dimensions; time dimension
%should be first, but script will fix if it's not
%signal_times - corresponding times to signal
%event_times - times of time-locked events (e.g. pips)
%before & after - # of samples before and after the event to calculate the
%phase consistency [must be positive numbers]
%freq_range - low and high cutoffs for bandpass filter
%fs - sampling frequency of signal
%order - size of the filter (optional; default: fs [1 second])
%
%OUT:
%phase_reset_strength - phase consistency at time points before and after event
%phases - event-triggered phases that were used to calculate consistency
%raw - event-triggered raw signal
%
%Joseph M. Stujenske, 2021
%
if nargin<7
    error('Requires at least 7 inputs.');
end
if nargin<8 || isempty(order)
    order=fs;
end
before=abs(before);
s_s=size(signal);
time_dim=s_s==length(signal_times);
if ~any(time_dim)
    error('Signal_times and signal do not align');
end
time_dim=find(time_dim,1,'first');
signal=permute(signal,[time_dim 1:time_dim-1 time_dim+1:length(s_s)]);
%get event-triggered signals, padded by 4000 samples due to filter
%inaccuracy on the edges
raw=event_triggered_signal(signal,signal_times,event_times,before+2000,after+2000);

%move dimensions around
n_events=length(event_times);
rep_dim=size(raw);
raw=permute(raw,[2:length(rep_dim) 1]);
rep_dim=length(rep_dim);
if n_events==1
    rep_dim=rep_dim+1;
end
%
%an inelegant solution to NaNs:
toremove=isnan(raw);
raw(toremove)=0;
[~,~,~,phases]=BandFilt_Order(raw,fs,order,freq_range(1),freq_range(2));
phases(toremove)=NaN;
phases=phases(2001:end-2000,:,:,:);
raw=raw(2001:end-2000,:,:,:);
phase_reset_strength=circ_r(phases,[],[],rep_dim);