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
toremove=ismember(gwfparams.spikeClusters,cids(cgs==0));
gwfparams.spikeTimes(toremove)=[];
gwfparams.spikeClusters(toremove)=[];
cids(cgs==0)=[];
cgs(cgs==0)=[];
wf = getWaveForms(gwfparams);
baseline_s=-100;
baseline=find((gwfparams.wfWin(1):gwfparams.wfWin(2))==baseline_s);
endplot_s=200;
endplot=find((gwfparams.wfWin(1):gwfparams.wfWin(2))==endplot_s);
wf.waveForms2=wf.waveForms;
M=memmapfile([day_dir,'\data_binary.bin']);
M=memmapfile([day_dir,'\data_binary.bin'],'Format',{'int16',[28 length(M.Data)/28/2],'data'});
A=M.Data.data(:,1:100000);
ranges=max(A,[],2)-min(A,[],2);
% noise_chan=ranges>nanmedian(double(ranges))*2;
noise_chan=ranges>min(double(ranges))*4;
wf.waveForms2(:,:,noise_chan,:)=NaN;
average_form=zeros(size(wf.waveForms2,[1 3 4]),'like',wf.waveForms2);
for a=1:size(wf.waveForms2,1)
average_form(a,:,:)=squeeze(nanmean(medfilt1(wf.waveForms2(a,:,:,:),5,[],4),2));
end
if plot_figures
figure;
    for n_elec=1:size(average_form,2)
        plot(reshape(squeeze(cat(3,average_form((cgs==2),n_elec,36:end),NaN*ones(sum(cgs==2),1,3)))',1,[])-(n_elec-1)*200,'k');hold on;
    end
end
size_wf=squeeze(min(nanmean(wf.waveForms2(:,:,:,baseline:endplot),2),[],4));
size_wf2=squeeze(max(abs(nanmean(wf.waveForms2(:,:,:,baseline:endplot),2)),[],4));
size_wf3=squeeze(nanmean(rms(wf.waveForms2(:,:,:,baseline:endplot),4),2));
% n_toplot=endplot-baseline+1;
if plot_figures;figure;end
spike_shape=zeros(size(average_form,1),endplot-baseline+1);
noise_unit=false(1,size(average_form,1));
for n=1:size(average_form,1)
in=find(size_wf(n,:)<-50);
if length(in)>14
    noise_unit(n)=true;
end
    [~,in]=sort(size_wf(n,:)./size_wf3(n,:),2,'ascend');
[~,in2]=sort(size_wf2(n,:),2,'ascend');
in2=in2(1:14);
in=in(1);
n_electoplot=min(length(in),1);
try
    baseLFP=average_form(n,in2,baseline+1:end);
    [~,variability]=sort(rms(baseLFP-nanmean(baseLFP,3),3),2,'ascend');
    
    backgroundLFP=squeeze(nanmedian(average_form(n,in2(variability(1:7)),baseline:end),2)); 
    temp=squeeze(nanmean(average_form(n,in(1:n_electoplot),baseline:end),2));
    lin_b=[backgroundLFP([1:60 160:end]) ones(203,1)]\temp([1:60 160:end]);
% spike_shape(n,:)=squeeze(nanmean(average_form(n,in(1:n_electoplot),baseline+1:endplot)-average_form(n,in(1:n_electoplot),start_t),2));
spike_shape(n,:)=temp(1:endplot-baseline+1)-(backgroundLFP(1:endplot-baseline+1)*lin_b(1)+lin_b(2));
catch
    keyboard
end
end
find_half=abs(spike_shape-spike_shape(:,100-baseline+1)/2);
[~,first_half]=min(find_half(:,1:100-baseline),[],2);
[~,second_half]=min(find_half(:,100-baseline+2:end),[],2);
halfwidth=(second_half+(100-baseline+1)-first_half)/3e1;
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
for nf=1:3
%     try
    spikeshape=data.clustermetrics.(fieldnames{nf}).SpikeShape;
    if ~isempty(spikeshape)
    spike_zero=spikeshape-nanmedian(spikeshape(:,1:50),2);
     data.clustermetrics.(fieldnames{nf}).SpikeShape_zero=spike_zero;
data.clustermetrics.(fieldnames{nf}).halfwidth2=[];
         for n_spike=1:size(spike_zero,1)
             try
            spike_trace=spike_zero(n_spike,:);
            spike_trace=BandFilt_Order(spike_trace,3e4,50,100,3000)';
             [maxval,in_temp]=max(spike_trace((spike_time-17):(spike_time+10)));
             [minval,in_temp2]=min(spike_trace((spike_time-17):(spike_time+10)));
             if (abs(maxval)>50 && abs(minval)>50 && in_temp<in_temp2 && maxval>abs(minval)*.5) || (maxval/abs(minval)>4)
                 spike_trace=-spike_trace;
             end
                 endsec=spike_trace(150:end);
                 linfit_b=[(1:length(endsec))'+149 ones(length(endsec),1)]\endsec';
                 spike_trace(120:end)=spike_trace(120:end)-((1:length(endsec)+30)-1)*linfit_b(1);
                 [spike_volt,spike_in]=min(spike_trace(spike_time-2:spike_time+15));
                 spike_in=spike_in+spike_time-3;
    av_afterspike=nanmedian(spike_trace(150:end));
                 fix_ratio=abs(spike_volt)/(av_afterspike-spike_volt);
                 spike_trace(spike_in:end)=(spike_trace(spike_in:end)-spike_volt)*fix_ratio+spike_volt;
                 spike_zero(n_spike,:)=spike_trace;
             spike_trace=resample(spike_trace,9e4,3e4);
%              spike_nt=length(spike_trace);
             [~,spike_bot]=nanmin(spike_trace((spike_time-2)*3:(spike_time+15)*3));
             spike_bot=spike_bot+(spike_time-2)*3-1;
    [~,first_half]=nanmin(abs(spike_trace(spike_time*3-30+1:spike_bot)-spike_trace(spike_bot)/2),[],2);
[~,second_half]=nanmin(abs(spike_trace(spike_bot+1:spike_bot+90)-spike_trace(spike_bot)/2),[],2);
halfwidth=(second_half+(spike_bot)-(first_half+spike_time*3-30))/9e1;
     data.clustermetrics.(fieldnames{nf}).halfwidth2(n_spike)=halfwidth(:)';
     spike_trace(isnan(spike_trace))=0;
lp=LowFilt_Order(spike_trace',9e4,100,1000)';
differential=diff([0 lp]);
peaks=find([false diff(differential>0)==-1]);
peaks(peaks>spike_bot+150)=[];
valleys=find([false diff(differential>0)==1]);
peak=peaks(find(peaks>spike_bot,1,'first'));
if isempty(peak)
    [~,peak]=min(abs(differential(spike_bot+second_half:spike_bot+150)));
    peak=peak+spike_bot+second_half-1;
    AUC=0;
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
data.clustermetrics.(fieldnames{nf}).troughpeaklatency(n_spike)=(peak-spike_bot)/9e1;
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