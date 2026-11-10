function build_ds23_data(dataRoot,workRoot,workers,cleanReferenceFile)
%BUILD_DS23_DATA 从本次自主跟踪结果和原始 IQ 计算五种统计量及绘图数据。
% 图的名称、色阶及标注单独定义；不加载冻结的检测结果或 I/Q 热图。
[cl,arrowTime,arrowDelay]=ds23_plot_style();
prn=[13 23 30 32;7 13 19 30];
contract=ds23_iq_contract();
rec=string({contract.receiver_folder});
meta=cell(1,3); setting=meta; tracks=meta;
if nargin<4,cleanReferenceFile='';end
useCleanFile=~isempty(char(cleanReferenceFile)); acq=meta; cleanShift=[74922938 55083021];
for k=1:3
 meta{k}=contract(k);
 path=fullfile(dataRoot,meta{k}.file);
 R=load(fullfile(workRoot,rec(k),'tracking.mat'),'tracks','settings');
 setting{k}=R.settings;setting{k}.fileName=char(path);
 tracks{k}=R.tracks;
 if useCleanFile
  acq{k}=load(fullfile(workRoot,rec(k),'acquisition.mat'),'acqResults');
 end
end
clear R
calPrn=unique(prn(:)).';
clean=cell(size(calPrn));
parfor (j=1:numel(calPrn),workers)
 tr=one_track(tracks{1},calPrn(j));
 clean{j}=derive_ds23_star(tr,setting{1},meta{1},[0 20],false,false);
end
caseData=cell(1,8);
parfor (j=1:8,workers)
 k=ceil(j/4);p=prn(k,mod(j-1,4)+1);
 tr=one_track(tracks{k+1},p);
 caseData{j}=derive_ds23_star(tr,setting{k+1},meta{k+1},[89.8 240],true,true);
 fprintf('Finished %s PRN %02d\n',rec(k+1),p);
end
out=fullfile(workRoot,'report');if ~isfolder(out),mkdir(out);end
D=struct('id',{},'title',{},'capture',{},'time',{},'value',{},'threshold',{}, ...
 'label',{},'thresholdLabel',{},'arrows',{});
B=struct([]);IQ=struct([]);results=struct([]);
methods={'sun','shang','deng','feng'};
titles={'Sun Q-SQM','Shang | Shape residual','Deng | Four-tap Q energy','Feng | Range curvature'};
labels={'5 s mean Q-SQM','Shape residual','Four-tap Q energy','Range curvature'};
gateLabels={'Clean 99th percentile','Nominal chi-square (160)','Nominal chi-square (4000)','Nominal chi-square (40)'};
for k=1:2
 name=char(extractBefore(rec(k+1),'_'));
 s=setting{k+1};M=meta{k+1};
 [ft,fz,fg,rho,accel,varF,nF]=clock_score(fullfile(workRoot,rec(k+1)),prn(k,:));
 [bt,bz,bg,br,ba,fit,detail]=adc_score(s.fileName,M);
 use=bt>=95&bt<=240;
 B(k).record=name; B(k).time=bt(use); B(k).value=bz(use); B(k).threshold=bg;
 B(k).rms=br(use); B(k).rmsMatch=ba(use); B(k).fit=fit;
 writetable(detail,fullfile(out,[name '_front_end.csv']));
 IQ(k).name=name;IQ(k).time=85:244;IQ(k).rms=nan(numel(IQ(k).time),2);
 n=round(.02*s.samplingFreq);
 for j=1:numel(IQ(k).time)
  first=round(IQ(k).time(j)*s.samplingFreq)-M.source_start_sample;
  x=read_iq(s.fileName,first,n);
  IQ(k).rms(j,:)=[sqrt(mean(real(x).^2)),sqrt(mean(imag(x).^2))];
 end
 for j=1:4
  i=(k-1)*4+j;p=prn(k,j);A=caseData{i};C=clean{calPrn==p};
  [sv,sg]=sun_qsqm(A.sunQ,A.promptI,A.valid,C.sunQ,C.promptI,C.valid,C.time>=5&C.time<20,5000,.01);
  shVariance=mean(A.cal)/160;sh=A.energy/shVariance;
  shGate=2*gammaincinv(1e-4,80,'upper');
  ix=A.time>=90;tap=A.fourQ(ix,:);v=A.valid(ix);v(1)=false;
  assert(all(v(2:4000))&&all(isfinite(tap(1:4000,:)),'all'),'Incomplete Deng baseline');
  [dv,dg]=deng_qcomposite(tap,v,tap(1:4000,:),4000,.01);
  stamps={A.time,A.time,A.time(ix),ft};scores={sv,sh,dv,fz(:,j)};gates=[sg shGate dg fg];
  C=A.capture;assert(isequal(C.t,(95:240).'));
  % 仅 DS2 PRN 13 有身份标注。其坐标不输入任何检测器。
  annotate=k==1&&p==13;
  if annotate&&useCleanFile
   prior=acq{1}.acqResults.carrFreq(p);
   if prior==0,prior=acq{k+1}.acqResults.carrFreq(p);end
   times=C.t;times(~isfinite(C.windowEnd))=NaN;
   C.clean=derive_clean_annotation(one_track(tracks{k+1},p),s,cleanReferenceFile,times,M.source_start_sample,cleanShift(k),prior);
  elseif annotate
   C.clean(ismember(C.t,arrowTime))=arrowDelay;
  end
  C.clean(~isfinite(C.windowEnd))=NaN;
  d=C.clean(ismember(C.t,arrowTime));
  annotate=annotate&&all(isfinite(d)&abs(d)<=2)&&d(end)<d(1);
  for m=1:4
   id=sprintf('%s_%s_prn%02d',methods{m},name,p);
   tt=stamps{m};zz=scores{m};inc=tt>=95&tt<=240;
   [tPlot,vPlot]=display_points(tt(inc),zz(inc));
   title=sprintf('%s | %s PRN %d',titles{m},upper(name),p);
   D(end+1)=struct('id',id,'title',title,'capture',C,'time',tPlot,'value',vPlot, ...
    'threshold',gates(m),'label',labels{m},'thresholdLabel',gateLabels{m},'arrows',m==1&&annotate); %#ok<AGROW>
  end
  if j==1,B(k).capture=repmat(C,1,4);end
  B(k).capture(j)=C;
  results(i).record=name;results(i).prn=p;
  results(i).sun_threshold=sg;results(i).shang_threshold=shGate;
  results(i).deng_threshold=dg;results(i).feng_threshold=fg;results(i).shang_variance=shVariance;
  results(i).deng_variance=var(sum(tap(1:4000,:),2),1);results(i).feng_variance=varF(j);
  results(i).feng_training_terms=nF(j);
  if k==2&&p==7
   ds3SunTime=A.time;ds3Sun=sv;ds3SunThr=sg;
   % 报告中近零分母：以同一历元的早迟 Q 和中心 I 重算。
   q=vecnorm(A.sunQ,2,2)./abs(A.promptI);q(~A.valid|~isfinite(q))=NaN;
   ii=find(A.time>=115&A.time<=160);[~,h]=max(q(ii),[],'omitnan');h=ii(h);
   note=table(A.time(h),A.promptI(h),norm(A.sunQ(h,:)),q(h),q(h)/5000, ...
    'VariableNames',{'time_s','prompt_i_adc','early_late_q_adc','epoch_score','contribution_to_mean'});
   writetable(note,fullfile(out,'prn7_high_plateau.csv'));
  end
 end
 rangeLabels=["time_s",compose("range_m_prn%02d",prn(k,:)), ...
  compose("range_accel_mps2_prn%02d",prn(k,:)),compose("feng_score_prn%02d",prn(k,:))];
 tableOut=array2table([ft rho accel fz],'VariableNames',cellstr(rangeLabels));
 writetable(tableOut,fullfile(out,[name '_pseudorange.csv']));
end
Z=scale_views(D);
file=fullfile(out,'figure_data.mat');save(file,'D','B','IQ','cl','Z','-v7');
writetable(struct2table(results),fullfile(out,'case_calibration.csv'));
render_figures(file,fullfile(out,'figures'));
diagnose_ds23_prn7(tracks{3},setting{3},meta{3},ds3SunTime,ds3Sun,ds3SunThr,out);
diagnose_ds23_residuals(tracks,setting,meta,out);
end

function tr=one_track(tracks,p)
%ONE_TRACK 使用本次捕获成功的指定通道；缺星时显式停止。
tr=tracks([tracks.PRN]==p);assert(isscalar(tr),'PRN %d missing or duplicated',p);
end

function Z=scale_views(D)
%SCALE_VIEWS 主图聚焦常态，另外画原始幅度。所有越出主图的点都在对应详图范围内。
Z=struct('id',{},'upper',{},'start',{},'stop',{},'peak',{});
for k=1:numel(D)
 E=D(k);v=E.value(isfinite(E.value));if isempty(v),continue;end
 ordinary=max(quantile(v,.90),E.threshold); peak=max(v);
 if strcmp(E.id,'sun_ds3_prn07')
  upper=1.5;
 elseif strcmp(E.id,'feng_ds3_prn07')
  upper=2000;
 elseif peak>4*ordinary
  upper=ceil(1.25*ordinary/10^floor(log10(ordinary)))*10^floor(log10(ordinary));
 else
  continue;
 end
 extreme=E.time(E.value>upper);assert(~isempty(extreme)&&upper>E.threshold);
 start=max(95,floor((min(extreme)-3)/5)*5);stop=min(240,ceil((max(extreme)+3)/5)*5);
 Z(end+1)=struct('id',E.id,'upper',upper,'start',start,'stop',stop,'peak',peak); %#ok<AGROW>
end
end

function [t,v]=display_points(t,v)
%DISPLAY_POINTS 高频数据保留40ms内首末点、极值、空段两侧。
t=t(:);v=v(:);if numel(t)<=8000,return;end
edges=[1;find(diff(floor(t/.04))~=0)+1;numel(t)+1];keep=false(size(t));
for j=1:numel(edges)-1
 i=edges(j):edges(j+1)-1;z=v(i);
 [~,lo]=min(z,[],'omitnan');[~,hi]=max(z,[],'omitnan');
 keep(i([1 end lo hi]))=true;
 b=find(diff(isfinite(z))~=0);keep(i([b;b+1]))=true;
end
t=t(keep);v=v(keep);
end

