function run_ds23_receiver(iq,meta,out,project,rx,isClean,acquireOnly)
%RUN_DS23_RECEIVER 固定 CU-SDR 捕获、跟踪、NAV，保存真实历元和导航观测。
% 调用者先核验 IQ 的大小和哈希。可选 acquireOnly 只运行实际捕获。
if nargin<7,acquireOnly=false;end
if ~isfolder(out),mkdir(out);end
s=initSettings();s.fileName=char(iq);s.dataType='int16';s.fileType=2;
s.samplingFreq=meta.sample_rate_hz;s.IF=meta.if_hz;
s.skipNumberOfBytes=0;s.resamplingflag=0;s.plotAcquisition=0;s.plotTracking=0;s.plotNavigation=0;
s.acqSatelliteList=1:32;
acqResults=acquisition(read_iq(iq,0,round(.05*s.samplingFreq)),s);
found=find(acqResults.carrFreq~=0);assert(~isempty(found),'No satellite acquired');
settings=s;save(fullfile(out,'acquisition.mat'),'acqResults','settings','-v7');
if acquireOnly,return;end
if isClean,duration=41.9;else,duration=160;end
settings.msToProcess=floor(duration*1000)-2;
settings.numberOfChannels=numel(found);
settings.navValidation='word_parity_and_issue_v1';
if isClean,settings.navValidation='native';end
channel=preRun(acqResults,settings);
batch=prepare_tracking(project,rx);addpath(batch,'-begin');clear tracking_logged
fid=fopen(iq,'rb','ieee-le');assert(fid>0); guard=onCleanup(@() fclose(fid));
[tracks,channel]=tracking_logged(fid,channel,settings);clear guard
save(fullfile(out,'tracking.mat'),'tracks','settings','channel','-v7.3');
if isClean
 [nav,eph]=postNavigation(tracks,settings);navAudit=[];
else
 batch=prepare_navigation(project,rx);addpath(batch,'-begin');clear postNavigation_checked decode_nav_subframe_checked
 [nav,eph,navAudit]=postNavigation_checked(tracks,settings);
end
save(fullfile(out,'navigation.mat'),'nav','eph','settings','navAudit','-v7');
if ~isempty(nav)&&isfield(nav,'currMeasSample')
 time_s=(nav.currMeasSample(:)+meta.source_start_sample)/settings.samplingFreq;
 writetable(table(time_s),fullfile(out,'positions.csv'));
end
if ~isClean,assert(isfile(fullfile(out,'positions.csv'))&&~isempty(nav),'NAV produced no measurement epochs');end
report=struct('source_iq_sha256',meta.sha256,'source_start_sample',meta.source_start_sample, ...
 'acquired_prns',found,'receiver_commit','0d506c7b3907113f9e371d9ff87be1aa83772917');
save(fullfile(out,'report.mat'),'report','-v7');
end
