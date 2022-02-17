function out=event_triggered_signal(signal,signal_times,event_times,before,after,average_dim,option)
%out=event_triggered_signal(signal,signal_times,event_times,before,after,average_dim,option)
%
%output a peri-event time matrix with the signal surround each of a series
%of events of interest
%you can get a peri-event time histogram by summing/averaging over the last
%dimension of the output matrix
%
%input:
%signal: n-dimensional signal (e.g. calcium-trace)
%signal_times: time stamps corresponding to signal
%event_times: time stamps corresponding to events
%before: number of frames/times before event to include
%after: number of frames/times after event to include
%average_dim: dimension to treat as time; default will be first dimension
%that is the same length as signal_times
%option: 'frame' (default) or 'time' -- uints for before and after
%'time' assumes that signal_times are monotonic increasing
%
%output:
%out: same dimensions as signal, except time dimension is the length of
%before+after+1, and a new first dimension is the length of event_times
%
%Joseph M. Stujenske, 2021
%
signal_times=double(signal_times(:)');
event_times=double(event_times(:)');
n_t=length(signal_times);
dim_size=(size(signal));
if nargin <6 || isempty(average_dim)
    average_dim=find(dim_size==n_t,1,'first');
end
if nargin <7 || isempty(option)
    option='frame';
end
if strcmp(option,'time')
    f_length=median(diff(signal_times));
    before=round(before/f_length);
    after=round(after/f_length);
end
signal=permute(double(signal),[average_dim,setdiff(1:length(dim_size),average_dim)]);
signal=reshape(signal,size(signal,1),[]);
dim_size2=size(signal);
% signal=signal(:)';
if any(event_times>signal_times(end))
    n=sum(event_times>signal_times(end));
    event_times(event_times>signal_times(end))=[];
    disp(['Ignoring ',num2str(n),' out of range events.'])
end
out_size=dim_size2;
out_size(1)=before+after+1;
out=zeros([length(event_times),out_size]);
for nCS=1:length(event_times)
    [~,event_frame]=min(abs(signal_times-event_times(nCS)));
    padlength_start=[];
    padlength_end=[];
    if event_frame-before<1
        padlength_start=1-(event_frame-before);
    end
    if event_frame+after>=length(signal)
        padlength_end=(event_frame+after)-length(signal);
    end
    temp=[NaN*ones(padlength_start,dim_size2(2));signal(max(event_frame-before,1):min(event_frame+after,length(signal)),:);NaN*ones(padlength_end,dim_size2(2))];

    out(nCS,:,:)=temp;
end
out_size=dim_size;
n_dim=length(out_size)+1;
out_size(average_dim)=before+after+1;
% out=permute(reshape(out,[length(event_times),dim_size2]),[2:n_dim 1]);
out=permute(out,[1 3:average_dim+1 2 average_dim+2:n_dim]);