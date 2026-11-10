function [score,gamma,ready,alarm]=bhuiyan_if(O,ref,p,cleanScore,mode)
%BHUIYAN_IF 原始ADC的Pearson统计量，O每列一个块，ref为冻结参考频数。
% cleanScore为对应分箱的有效基线统计量；gaussian用修正的上尾，
% empirical用同一分箱的经验分位。返回线性值，绘图时可转dB。
assert(size(O,1)==numel(ref)&&all(isfinite(O)&O>=0,'all')&&all(isfinite(ref(:))&ref(:)>=0)&&sum(ref)>0&&p>0&&p<1);
P=ref(:)/sum(ref);use=P>0;E=P*sum(O,1);x=O(use,:);e=E(use,:);
ready=sum(O,1)>0&~any(O(~use,:)>0,1);score=nan(1,size(O,2));
score(ready)=sum((x(:,ready)-e(:,ready)).^2./e(:,ready),1);
h=cleanScore(isfinite(cleanScore));gamma=NaN;
if ~isempty(h)
 switch mode
  case 'gaussian',gamma=mean(h)+sqrt(2)*std(h,1)*erfcinv(2*p);
  case 'empirical',gamma=quantile(h,1-p);
  otherwise,error('Unknown threshold mode');
 end
end
alarm=ready&score>=gamma;
end
