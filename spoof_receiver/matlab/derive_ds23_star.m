function A=derive_ds23_star(tr,s,M,limits,fullTaps,makeMap)
%DERIVE_DS23_STAR 逐历元重相关；完整性检查只看实际样本边界、NCO 和连续性。
% limits 为原录制时刻上下界；fullTaps 为 Shang / Deng 生成整峰能量及四点 Q。
% makeMap 另按整秒构造 20 个连续积分的 I/Q 功率图。
a=double(tr.absoluteSample(:));b=double(tr.endSample(:));fs=s.samplingFreq;
t=(b+M.source_start_sample)/fs;n=b-a;
v=isfinite(a)&isfinite(b)&a>=0&b<=M.sample_count&n>0&isfinite(tr.codeFreq(:)) & ...
 isfinite(tr.carrFreq(:))&isfinite(tr.remCodePhase(:))&isfinite(tr.remCarrPhase(:));
% 首个实际完整积分没有前驱，不因缺少前驱而丢掉首个 5 秒窗口。
continuous=[true;a(2:end)==b(1:end-1)&abs(diff(t)-.001)<.0001];
idx=find(t>=limits(1)&t<=limits(2)); m=numel(idx); code=generateCAcode(tr.PRN);
q=nan(m,2);ip=nan(m,1);energy=ip;four=nan(m,4);cal=nan(90,1);
for base=1:256:m
 j=base:min(base+255,m);sel=idx(j);g=v(sel);if ~any(g),continue;end
 first=min(a(sel(g))); last=max(b(sel(g)));
 raw=read_iq(s.fileName,first,last-first); U=nan(81,numel(j));
 for u=find(g(:)).'
  ep=sel(u);x=raw(a(ep)-first+(1:n(ep)));
  z=recompute_correlators(s.fileName,tr,ep,s,[-.4 0 .4],x,code)/n(ep);
  q(j(u),:)=imag(z([1 3]));ip(j(u))=real(z(2));
  if fullTaps
   w=recompute_correlators_folded(s.fileName,tr,ep,s,x)/n(ep);
   U(:,u)=w.';four(j(u),:)=imag(w([31 36 46 51]));
  end
 end
 if fullTaps,energy(j(g))=shang_multicorrelator(U(:,g),20,1,1e-4).';end
end
if fullTaps
 grid=86:.1:94.9; U=nan(81,numel(grid));
 for j=1:numel(grid)
  k=find(t<=grid(j),1,'last');assert(~isempty(k)&&v(k),'Missing Shang baseline frame');
  U(:,j)=recompute_correlators_folded(s.fileName,tr,k,s).'/n(k);
 end
 cal=shang_multicorrelator(U,20,1,1e-4).';
end
A=struct('time',t(idx),'sunQ',q,'promptI',ip,'valid',v(idx)&continuous(idx), ...
 'energy',energy,'fourQ',four,'cal',cal,'capture',[]);
if makeMap
 times=(max(95,ceil(limits(1))):min(240,floor(limits(2)))).';d=(-40:40)/20;
 I=nan(numel(times),81,'single');Q=I;peak=nan(size(times));goodEnd=peak;
 for j=1:numel(times)
  target=round(times(j)*fs)-M.source_start_sample;
  k=find(b<=target & v,1,'last');
  if isempty(k)||k<20,continue;end
  ids=k-19:k;
  if target-b(k)>1.1*fs/1000||~all(v(ids))||~all(continuous(ids(2:end))),continue;end
  [i,q,peak(j)]=profile(s.fileName,tr,s,ids,code,d,false);
  I(j,:)=single(10*log10(i));Q(j,:)=single(10*log10(q));goodEnd(j)=t(k);
 end
 A.capture=struct('t',times,'d',d(:),'I',I,'Q',Q,'clean',nan(size(times)), ...
     'peak',peak,'prn',tr.PRN,'windowEnd',goodEnd);
end
end

function [I,Q,peak,err,trerr]=profile(file,tr,s,ix,code,d,verify)
%PROFILE 原IQ连续读一次；20个1ms历元各自相干，分支功率再做非相干平均。
st=double(tr.absoluteSample(ix)); en=double(tr.endSample(ix)); count=en-st;
x=read_iq(file,st(1),en(end)-st(1)); U=zeros(numel(ix),numel(d)); P0=zeros(size(ix));
err=0;trerr=0;
for j=1:numel(ix)
    n=count(j); b=x(st(j)-st(1)+(1:n)); ep=ix(j);
    % 增大物理延迟等价于减小本地码相位；所有点共用同一跟踪载波。
    z=taps(file,tr,ep,s,b,code,-d); U(j,:)=z/n;
    P0(j)=sum(abs(b).^2)/n^2;
    if verify && (j==1||j==numel(ix))
        v=recompute_correlators(file,tr,ep,s,-d,b,code);
        err=max(err,max(abs(z-v))/max(1,max(abs(v))));
        p=complex(tr.I_P(ep),tr.Q_P(ep));
        trerr=max(trerr,abs(p-z(d==0))/max(1,abs(p)));
    end
end
bg=mean(P0); assert(bg>0);
I=mean(real(U).^2,1)/bg; Q=mean(imag(U).^2,1)/bg;
[~,i]=max(I+Q); peak=d(i);
end

function z=taps(file,tr,ep,s,raw,code,offset)
%TAPS 按本地C/A码的码片边界，对每个延迟的样本前缀和分段求和。
% 每片贡献仍是码符号乘本片内全部原始样本；边界数值不确定时用逐样本公式。
n=numel(raw); f=s.samplingFreq; step=tr.codeFreq(ep)/f;
base=raw.*exp(-1i*(2*pi*tr.carrFreq(ep)*((0:n-1)/f)+tr.remCarrPhase(ep)));
a=tr.remCodePhase(ep)+offset;
b=(floor(min(a))-1:ceil(max(a)+(n-1)*step)).';
edge=(b-a)/step; cut=min(max(floor(edge)+1,0),n);
prefix=[0 cumsum(base)]; sums=diff(prefix(cut+1),1,1);
z=code(mod(b(2:end)-1,s.codeLength)+1)*sums;
near=any(edge>=-1e-7 & edge<=n-1+1e-7 & abs(edge-round(edge))<1e-7,1);
if any(near),z(near)=recompute_correlators(file,tr,ep,s,offset(near),raw,code);end
end
