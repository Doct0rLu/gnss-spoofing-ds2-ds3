function [t,score,g,rho,a,s2,N]=clock_score(folder,prns)
%CLOCK_SCORE 用固定采样钟代替位置解算反馈的钟修正，保留电文不可用的时刻。
S=load(fullfile(folder,'navigation.mat'),'nav','settings'); V=S.nav; s=S.settings;
T=readtable(fullfile(folder,'positions.csv')); t=T.time_s(:);
assert(numel(t)==size(V.rawP,2)); x=V.currMeasSample(:); rho=nan(numel(t),numel(prns));
for j=1:numel(prns)
 tx=nan(size(t)); raw=tx;
 for k=1:numel(t)
  c=find(V.PRN(:,k)==prns(j),1);
  if ~isempty(c),tx(k)=V.transmitTime(c,k);raw(k)=V.rawP(c,k);end
 end
 k=find(isfinite(raw)&isfinite(tx),1);if isempty(k),continue;end
 rho(:,j)=raw(k)+s.c*((x-x(k))/s.samplingFreq-(tx-tx(k)));
end
dt=median(diff(t)); assert(abs(dt-.5)<1e-6);
a=nan(size(rho));a(3:end,:)=diff(rho,2,1)/dt^2;
step=abs(diff(t)-dt)<1e-6;a(~[false;false;step(1:end-1)&step(2:end)],:)=NaN;
train=t>=91&t<99; tri=[false;false;train(3:end)&train(2:end-1)&train(1:end-2)];
N=sum(isfinite(a(tri,:)),1);s2=var(a(tri,:),1,1,'omitnan');s2(N<3)=NaN;
[score,g,~,~,~,check]=feng_cd2p(t,rho,zeros(size(rho)),s2,40,1e-5);
assert(max(abs(a-check),[],'all','omitnan')<1e-8);
end

