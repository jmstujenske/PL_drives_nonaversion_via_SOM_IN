function data=addwaveforms(data,plot_figures)
if nargin<2 || isempty(plot_figures)
    plot_figures=false;
end
day_dir=data.rootdir;
gwfparams.dataDir=day_dir;
gwfparams.fileName='data_binary.bin';
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 28;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-100 201];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =readNPY([gwfparams.dataDir,'\spike_times.npy']);

try
[cids,cgs]=getclustermarkings(day_dir);
cl=readNPY([day_dir,'\spike_clusters.npy']);
catch
    cl=[];
end
gwfparams.spikeClusters=cl;

%remove noise clusters
toremove=ismember(gwfparams.spikeClusters,cids(cgs==0));
gwfparams.spikeTimes(toremove)=[];
gwfparams.spikeClusters(toremove)=[];
cids(cgs==0)=[];
cgs(cgs==0)=[];

%get waveforms (using script form kilosort analysis scripts)
wf = getWaveForms(gwfparams);
baseline_s=-100;
baseline=find((gwfparams.wfWin(1):gwfparams.wfWin(2))==baseline_s);
endplot_s=200;
endplot=find((gwfparams.wfWin(1):gwfparams.wfWin(2))==endplot_s);

wf.waveForms2=wf.waveForms;

%Memory map first just to figure out how many datapoints there are; alternatively could use fread, fseek to eof, ftell to get size
M=memmapfile([day_dir,'\data_binary.bin']);

%length(M.data) will be twice as long as the number of values because the values take up two bytes each
%Now memory map with knowledge of the format
M=memmapfile([day_dir,'\data_binary.bin'],'Format',{'int16',[gwfparams.nCh length(M.Data)/gwfparams.nCh/2],'data'});

%Arbitrarily take the first 100,000 data points to identify noisy channels
A=M.Data.data(:,1:min(100000,size(M.Data.data,2)));

%Noisy channels have very large positive and negative values; let's say that noisy channels will have 4 times larger range of values than the smallest amplitude channel
ranges=max(A,[],2)-min(A,[],2);
% noise_chan=ranges>nanmedian(double(ranges))*2;
noise_chan=ranges>min(double(ranges))*4;
%

%Remove the noise channels from waveForms2
wf.waveForms2(:,:,noise_chan,:)=NaN;

%Preallocate matrix for the waveforms
average_form=zeros(size(wf.waveForms2,[1 3 4]),'like',wf.waveForms2);

%Get the average in every channel across spikes for each cell
for a=1:size(wf.waveForms2,1)
average_form(a,:,:)=squeeze(nanmean(medfilt1(wf.waveForms2(a,:,:,:),5,[],4),2));
end

%Plot the average waveform across channels for each "good" cell
if plot_figures
figure;
    for n_elec=1:size(average_form,2)
        plot(reshape(squeeze(cat(3,average_form((cgs==2),n_elec,36:end),NaN*ones(sum(cgs==2),1,3)))',1,[])-(n_elec-1)*200,'k');hold on;
    end
end

%Different ways of assessing sizes of the waveforms
size_wf=squeeze(min(nanmean(wf.waveForms2(:,:,:,baseline:endplot),2),[],4));
size_wf2=squeeze(max(abs(nanmean(wf.waveForms2(:,:,:,baseline:endplot),2)),[],4));
size_wf3=squeeze(nanmean(rms(wf.waveForms2(:,:,:,baseline:endplot),4),2));
% n_toplot=endplot-baseline+1;

%Open new figure if plotting
if plot_figures;figure;end

%Preallocate matrix for spike shapes
spike_shape=zeros(size(average_form,1),endplot-baseline+1);

%Preallocate matrix for excluding units
noise_unit=false(1,size(average_form,1));

%Loop through cells
for n=1:size(average_form,1)

%A way that I found to empirally identify noise (large deflections across over half the channels)
%This may not be a generalizable solution
in=find(size_wf(n,:)<-50);
if length(in)>floor(gwfparams.nCh/2)
    noise_unit(n)=true;
end

%sort the channels by waveform size
[~,in]=sort(size_wf(n,:)./size_wf3(n,:),2,'ascend'); %this is a negative number so this orders from largest to smallest deflection
[~,in2]=sort(size_wf2(n,:),2,'ascend'); %these numbers are positive, so order is smallest to largest
in2=in2(1:ceil(gwfparams.nCh/2)); %pick the half of channels with the smallest deflection
in=in(1); %Channel with largest evoked change

%What follows below is a quick waveform extraction and calculation that doesn't work very well
%In this case, n will just be 1, but could change line above to use the top n channels
n_electoplot=min(length(in),1);
try

    %baseline LFP is inferred from half of channels that are least likely to contain the unit
    baseLFP=average_form(n,in2,baseline+1:end); %get baseline for unit n
    [~,variability]=sort(rms(baseLFP-nanmean(baseLFP,3),3),2,'ascend'); %find channels with smallest baseline LFP fluctuations
    
    backgroundLFP=squeeze(nanmedian(average_form(n,in2(variability(1:ceil(gwfparams.nCh/4))),baseline:end),2)); %take median of these channels to get the LFP that the spike is superimposed on
    temp=squeeze(nanmean(average_form(n,in(1:n_electoplot),baseline:end),2)); %now take the channels containing the spike
    lin_b=[backgroundLFP([1:60 160:end]) ones(203,1)]\temp([1:60 160:end]); %get a linear regressor to fit the LFP background to the channel with the spike, b/c amplitude not guaranteed to be the same btwn channels
% spike_shape(n,:)=squeeze(nanmean(average_form(n,in(1:n_electoplot),baseline+1:endplot)-average_form(n,in(1:n_electoplot),start_t),2));
spike_shape(n,:)=temp(1:endplot-baseline+1)-(backgroundLFP(1:endplot-baseline+1)*lin_b(1)+lin_b(2)); %remove the background LFP
catch
    keyboard
end
end

%find the half-way point of the spike on the left
find_half=abs(spike_shape-spike_shape(:,100-baseline+1)/2);
[~,first_half]=min(find_half(:,1:100-baseline),[],2);

%... and right
[~,second_half]=min(find_half(:,100-baseline+2:end),[],2);

%half-width should be the time difference between these points
halfwidth=(second_half+(100-baseline+1)-first_half)/3e1;

%Now assign units to the correct groupings
fieldnames={'noise','MUA','good'};
counts=[1 1 1];
allcl=unique(cl);
cl_count=1;
for a=1:3
data.clustermetrics.(fieldnames{a}).Halfwidth=[];
data.clustermetrics.(fieldnames{a}).SpikeShape=[];
end
for a=1:length(cids)
    if ismember(cids(a),allcl)
data.clustermetrics.(fieldnames{cgs((cids==cids(a)))+1}).Halfwidth(counts(cgs((cids==cids(a)))+1))=halfwidth(cl_count);
    data.clustermetrics.(fieldnames{cgs((cids==cids(a)))+1}).SpikeShape(counts(cgs((cids==cids(a)))+1),:)=spike_shape(cl_count,:);
    data.clustermetrics.(fieldnames{cgs((cids==cids(a)))+1}).Noise(counts(cgs((cids==cids(a)))+1))=noise_unit(cl_count);
    cl_count=cl_count+1;
    else
    data.clustermetrics.(fieldnames{cgs((cids==cids(a)))+1}).Halfwidth(counts(cgs((cids==cids(a)))+1))=NaN;
    data.clustermetrics.(fieldnames{cgs((cids==cids(a)))+1}).Noise(counts(cgs((cids==cids(a)))+1))=true;
    data.clustermetrics.(fieldnames{cgs((cids==cids(a)))+1}).SpikeShape(counts(cgs((cids==cids(a)))+1),:)=NaN;
    end
    counts(cgs(a)+1)=counts(cgs((cids==cids(a)))+1)+1;
end

%Now let's do this again but better.
for nf=1:3
%     try
    spikeshape=data.clustermetrics.(fieldnames{nf}).SpikeShape;
    if ~isempty(spikeshape)
    
    %Take the spike shape that was extracted previously and zero to before the spike
    spike_zero=spikeshape-nanmedian(spikeshape(:,1:50),2);
    data.clustermetrics.(fieldnames{nf}).SpikeShape_zero=spike_zero;
    
    %Pre-initialize a new half-width calculation that worsks better
    data.clustermetrics.(fieldnames{nf}).halfwidth2=[];
         for n_spike=1:size(spike_zero,1) %loop through units
             try
            spike_trace=spike_zero(n_spike,:); %get the trace for unit n_spike
            spike_trace=BandFilt_Order(spike_trace,3e4,50,100,3000)'; %bandpass filter to get rid of low frequency and very high frequency content
             [maxval,in_temp]=max(spike_trace((spike_time-17):(spike_time+10))); %find the maximum time 
             [minval,in_temp2]=min(spike_trace((spike_time-17):(spike_time+10))); % find the minimum time
             
             %Deal with the fact that some spikes are positive-going and flip them-- this MOSTLY identified these correctly
             %This needs some work because it sometimes gets it wrong
             if (abs(maxval)>50 && abs(minval)>50 && in_temp<in_temp2 && maxval>abs(minval)*.5) || (maxval/abs(minval)>4)
                 spike_trace=-spike_trace;
             end
             
             %Find the new baseline after the spike -- due to the underlying LFP this could be above or below zero
                 endsec=spike_trace(150:end);
                 
                 %Fit a line across the spike so that we can then make the baseline approximately flat
                 linfit_b=[(1:length(endsec))'+149 ones(length(endsec),1)]\endsec';
                 spike_trace(120:end)=spike_trace(120:end)-((1:length(endsec)+30)-1)*linfit_b(1);
                 
                 %Try to figure out when the spike occurred
                 [spike_volt,spike_in]=min(spike_trace(spike_time-2:spike_time+15));
                 spike_in=spike_in+spike_time-3; %deal with the fact that we did this on a sub-range
                 
                 %The baseline is usually still not flat, so now let's rescale things to make it approximately so
                 av_afterspike=nanmedian(spike_trace(150:end));
                 fix_ratio=abs(spike_volt)/(av_afterspike-spike_volt);
                 spike_trace(spike_in:end)=(spike_trace(spike_in:end)-spike_volt)*fix_ratio+spike_volt;
                 spike_zero(n_spike,:)=spike_trace;
                 
                 %upsample the spike trace
                spike_trace=resample(spike_trace,9e4,3e4);
%              spike_nt=length(spike_trace);
             
              %Try to figure out when the spike occurred
             [~,spike_bot]=nanmin(spike_trace((spike_time-2)*3:(spike_time+15)*3)); %re-estimate the spike time in the upsampled situation
             spike_bot=spike_bot+(spike_time-2)*3-1; %make it relative to the start because we did this on a small snippet
    [~,first_half]=nanmin(abs(spike_trace(spike_time*3-30+1:spike_bot)-spike_trace(spike_bot)/2),[],2); %find the left size of the half-width
[~,second_half]=nanmin(abs(spike_trace(spike_bot+1:spike_bot+90)-spike_trace(spike_bot)/2),[],2); %and the right size
halfwidth=(second_half+(spike_bot)-(first_half+spike_time*3-30))/9e1; %get the half-width as before
     data.clustermetrics.(fieldnames{nf}).halfwidth2(n_spike)=halfwidth(:)';
     spike_trace(isnan(spike_trace))=0;
     
     %1000 Hz low pass filter to try to get the after hyperpolarization
lp=LowFilt_Order(spike_trace',9e4,100,1000)';
differential=diff([0 lp]); %get the change in the low pass signal [the 0 is just to pad so the size stays the same]
peaks=find([false diff(differential>0)==-1]); %find positive inflection points
peaks(peaks>spike_bot+150)=[]; %ignore inflection points that occur too late
valleys=find([false diff(differential>0)==1]); %find negative inflection points
peak=peaks(find(peaks>spike_bot,1,'first'));  %find the first peak after the spike

%calculate area under the curve for after hyperpolarization and repolarization time (when it goes back to zero baseline
if isempty(peak)
    [~,peak]=min(abs(differential(spike_bot+second_half:spike_bot+150))); %there is no peak, so just find where things flatten the most
    peak=peak+spike_bot+second_half-1; %we calculated this on a snippet, so put relative to start
    AUC=0; %no after hyperpolarization, so it is zero
    second_half2=find(spike_trace(peak+1:end)>0 & spike_trace(peak:end-1)<0,1,'first')+peak;
    data.clustermetrics.(fieldnames{nf}).RepolarizationTime(n_spike)=(second_half2-peak)/9e1;
else
    afterbase=spike_trace(valleys(find(valleys>peak,1,'first')));
    second_half2=valleys(find(valleys>peak,1,'first'));
    first_half2=find(lp(spike_bot:peak)<=afterbase,1,'last');
    data.clustermetrics.(fieldnames{nf}).SpikeShape_zero(n_spike,:)=resample(lp,3e4,9e4);
    spike_norm=abs(nansum(lp(spike_time*3-30+first_half:second_half+(spike_bot))));
    spike_trace_rms=lp./(spike_norm);
    AUC=sum((spike_trace_rms(first_half2+spike_bot-1:second_half2)-afterbase/spike_norm)/9e1);
    data.clustermetrics.(fieldnames{nf}).RepolarizationTime(n_spike)=(second_half2-peak)/9e1;
end
data.clustermetrics.(fieldnames{nf}).areaunderpeak(n_spike)=AUC;
%another metric: spike to peak latency
data.clustermetrics.(fieldnames{nf}).troughpeaklatency(n_spike)=(peak-spike_bot)/9e1;

%size of the spike
data.clustermetrics.(fieldnames{nf}).amplitude(n_spike)=spike_trace(spike_bot);
             catch
         data.clustermetrics.(fieldnames{nf}).halfwidth2(n_spike)=NaN;
         data.clustermetrics.(fieldnames{nf}).areaunderpeak(n_spike)=NaN;
         data.clustermetrics.(fieldnames{nf}).troughpeaklatency(n_spike)=NaN;
         disp('no spikes.')
             end
         end
    else
                 data.clustermetrics.(fieldnames{nf}).halfwidth2=[];
         data.clustermetrics.(fieldnames{nf}).areaunderpeak=[];
         data.clustermetrics.(fieldnames{nf}).troughpeaklatency=[];
         data.clustermetrics.(fieldnames{nf}).SpikeShape_zero=[];
    end
end
keyboard
end
