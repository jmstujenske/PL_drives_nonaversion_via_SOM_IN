function data=process_post_phy(day_dir,animal_data,exclude_electrodes)
if nargin<3 || isempty(exclude_electrodes)
    exclude_electrodes=[];
end
try
[cids,cgs]=getclustermarkings(day_dir);
spt=readNPY([day_dir,'\spike_times.npy']);
cl=readNPY([day_dir,'\spike_clusters.npy']);
features=readNPY([day_dir,'\pc_features.npy']);
catch
    spt=[];
    cl=[];
end
data=struct();
data.rootdir=day_dir;
fieldnames={'noise','MUA','good'};
for a=max(cl)'+1;end

if ~isempty(spt)
data.units.good=cell(1,sum(cgs==2));
data.units.MUA=cell(1,sum(cgs==1));
data.units.noise=cell(1,sum(cgs==0));
counts=[1 1 1];
for a=1:length(cids)
    spike_times=spt(cl==cids(a));
    data.units.(fieldnames{cgs(a)+1}){counts(cgs(a)+1)}=spike_times;
    [IsolDis,Lratio]=IsolationDistance(double(reshape(features,size(features,1),[])),find(cl==cids(a)),[],find(cgs==0));
    data.clustermetrics.(fieldnames{cgs(a)+1}).IsolationDistance(counts(cgs(a)+1))=IsolDis;
    data.clustermetrics.(fieldnames{cgs(a)+1}).Lratio(counts(cgs(a)+1))=Lratio;
    counts(cgs(a)+1)=counts(cgs(a)+1)+1;
end
else
    data.units.good=[];
    data.units.MUA=[];
    data.units.noise=[];
end
try
data=addwaveforms(data);
end
try
electrodenumbers=vertcat(animal_data.NS5.ElectrodesInfo.ElectrodeID);
catch
    disp('Error. Entering debug mode.');
    keyboard
end
trunc_data=animal_data.NS5.Data(~ismember(electrodenumbers,exclude_electrodes),:);
trunc_data=nanmedian(trunc_data,1);
[trunc_data] = removeLineNoise_SpectrumEstimation(trunc_data, 30000, 'LF = 60 NH = 5');
trunc_data=lowpass(single(trunc_data),300,30000);
trunc_data=trunc_data(1:15:end);
% PL=animal_data.NS5.Data(electrodenumbers==32,:);
% 
% [PL] = removeLineNoise_SpectrumEstimation(PL, 30000, 'LF = 60 NH = 5');
% PL=lowpass(single(PL),300,30000);
% PL=PL(1:15:end);
% data.LFP_PL=int16(PL*4);
% mmf = memmapfile([day_dir,'\data_binary.bin'], 'Format', {'int16', size(animal_data.NS5.Data), 'x'});
% trunc_data=mmf.Data.x(:,1:15:end);
trunc_data=int16(trunc_data*4);%downsample to 2000 Hz
data.LFP_ts=1:15:size(animal_data.NS5.Data,2);
data.LFP_med=nanmedian(trunc_data);
data.LFP_full=trunc_data;
clear trunc_data;
data.Video_ts=animal_data.NEV.Data.VideoSync.TimeStamp;
LaserOnOfforSoundOff=animal_data.NEV.Data.SerialDigitalIO.UnparsedData==65534 | animal_data.NEV.Data.SerialDigitalIO.UnparsedData==65535;
LaserOnOfforSoundOff_ts=animal_data.NEV.Data.SerialDigitalIO.TimeStamp(LaserOnOfforSoundOff);
LaserOnOfforSoundOff_ev=find(LaserOnOfforSoundOff);
index=diff(LaserOnOfforSoundOff_ts)/30000;LaserOn=[index>30 & index<40 false];LaserOff=[false index>4 & index<6];
ToneOff=~(LaserOn | LaserOff);
LaserOn_ev=LaserOnOfforSoundOff_ev(LaserOn);
LaserOff_ev=LaserOnOfforSoundOff_ev(LaserOff);
ToneOff_ev=LaserOnOfforSoundOff_ev(ToneOff);
ToneOn_ev=setdiff(1:length(animal_data.NEV.Data.SerialDigitalIO.TimeStamp),sort([LaserOn_ev;LaserOff_ev;ToneOff_ev]));

data.ToneOn_ts=animal_data.NEV.Data.SerialDigitalIO.TimeStamp(ToneOn_ev);
data.ToneOff_ts=animal_data.NEV.Data.SerialDigitalIO.TimeStamp(ToneOff_ev);
data.LaserOn_ts=animal_data.NEV.Data.SerialDigitalIO.TimeStamp(LaserOn_ev);
data.LaserOff_ts=animal_data.NEV.Data.SerialDigitalIO.TimeStamp(LaserOff_ev);
data.shockfreq=animal_data.shock;
data.ToneswithLaser=logical(histc(data.LaserOn_ts,[data.ToneOn_ts-6*3e5 data.ToneOff_ts(end)]));
data.ToneswithLaser=data.ToneswithLaser(1:end-1);
TwoKHz=mod(animal_data.cueseq,2)>0;
data.ToneFreq=zeros(1,length(data.ToneOn_ts));
data.ToneFreq(TwoKHz)=2;
data.ToneFreq(~TwoKHz)=8;
data.CueSeq=animal_data.cueseq;
