function [t,z,g,r,rAdj,fit,detail]=adc_score(path,M)
%ADC_SCORE 用同一初始窗标定两段录制，另做幅度等化以定位直方图变化的来源。
fs=M.sample_rate_hz;n=round(.1*fs);start=M.source_start_sample/fs;
assert(strcmp(M.status,'complete')&&start==85);
base=zeros(10*n,1,'single');refms=zeros(5,1);
for s=0:4
 x=read_iq(path,s*fs,n);base((2*s*n+1):(2*s+1)*n)=real(x);base(((2*s+1)*n+1):(2*s+2)*n)=imag(x);
 refms(s+1)=mean(real(x).^2);
end
edges=double(unique(quantile(base,linspace(0,1,65))'));edges(1)=-inf;edges(end)=inf;
counts=histcounts(base,edges).';clear base x
assert(numel(counts)==64&&min(counts)*n/sum(counts)>5);
t=(85:244).'+.1;r=nan(size(t));raw=zeros(64,numel(t)); adjusted=raw;baselineRms=sqrt(mean(refms));
for k=1:numel(t)
 u=real(read_iq(path,(k-1)*fs,n));r(k)=sqrt(mean(u.^2));
 raw(:,k)=histcounts(u,edges).';
 adjusted(:,k)=histcounts(u*(baselineRms/r(k)),edges).';
end
[v,~,ready]=bhuiyan_if(raw,counts,.01,[],'empirical'); assert(all(ready)&&all(v>0));
fit=v(t<99);assert(numel(fit)==14);
[v,g0]=bhuiyan_if(raw,counts,.01,fit,'empirical');
w=bhuiyan_if(adjusted,counts,.01,fit,'empirical');
z=10*log10(v(:)); g=10*log10(g0); rAdj=10*log10(w(:));
detail=table(t,r,z,rAdj,'VariableNames',{'time_s','rms_adc','score_dB','score_dB_at_baseline_rms'});
end
