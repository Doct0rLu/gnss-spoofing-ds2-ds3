function clean=derive_clean_annotation(track,settings,cleanFile,times,sourceStart,offset,prior)
%DERIVE_CLEAN_ANNOTATION 仅重算报告标注所用的独立干净峰，不输入任何检测器。
% cleanFile 为从 0 s 开始的官方 cleanStatic 完整录制；offset 为配对复采样偏移。
% times 为攻击录制的原时间，prior 为捕获得到的载波频率初值（Hz）。
fs=settings.samplingFreq;N=round(fs/1000);step=1.023e6/fs;
ca=generateCAcode(track.PRN);code=ca(ceil((1:N)*step));clean=nan(size(times));
a=double(track.absoluteSample(:));b=double(track.endSample(:));
for k=1:numel(times)
 if ~isfinite(times(k)),continue;end
 pos=round((times(k)-.02)*fs);
 x=read_iq(cleanFile,pos+offset,20*N);
 coarse=prior+(-1000:250:1000);P=acquisition_power_map(x(1:5*N),code,fs,coarse);
 [~,p]=max(P(:));[i,~]=ind2sub(size(P),p);
 for width=[150 500]
  f=coarse(i)+(-width:25:width);P=acquisition_power_map(x,code,fs,f);
  [~,p]=max(P(:));[u,lag]=ind2sub(size(P),p);
  if u>1&&u<numel(f),break;end
 end
 first=pos-sourceStart;
 ep=find(isfinite(b)&a<=first&b>first,1,'last');
 if isempty(ep)||u==1||u==numel(f),continue;end
 d=a(ep)-track.remCodePhase(ep)/(track.codeFreq(ep)/fs)-(first+lag-2);
 clean(k)=-(mod(d+N/2,N)-N/2)*step;
end
end
