function [score,gamma,alarm,gaussian]=sun_qsqm(Q,IP,valid,cleanQ,cleanIP,cleanValid,calMask,w,p)
%SUN_QSQM 早迟Q归一化幅度；Q各列为±0.4 chip，IP为Prompt的I。
% clean输入与攻击段独立；calMask选取干净滑窗的结束历元。
% valid和cleanValid仅描述样本完整和连续；不按相位质量或载噪比筛选。
assert(w>=1&&w==fix(w)&&p>0&&p<1);
[score,~]=window(Q,IP,valid,w);
[h,q]=window(cleanQ,cleanIP,cleanValid,w);
cal=logical(calMask(:));base=h(cal&isfinite(h));raw=q(cal&isfinite(q));
gamma=NaN;gaussian=NaN;
if ~isempty(base),gamma=quantile(base,1-p);end
if ~isempty(raw)
 sigma=mean(raw)/sqrt(pi/2);
 gaussian=sigma*sqrt(pi/2)+sigma*sqrt((4-pi)/w)*erfcinv(2*p);
end
alarm=isfinite(score)&score>gamma;
end
function [s,q]=window(Q,IP,valid,w)
% 对连续有效历元求因果均值；无效位置中断窗口。
IP=IP(:);assert(isequal(size(Q),[numel(IP) 2]));
a=abs(IP);v=logical(valid(:)) & all(isfinite(Q),2) & isfinite(a);
v=v & a>0;
q=nan(size(a));q(v)=hypot(Q(v,1),Q(v,2))./a(v);
v=v & isfinite(q); q(~v)=NaN;
x=q;x(~v)=0;s=movsum(x,[w-1 0])/w;s(movsum(double(v),[w-1 0])~=w)=NaN;
end
