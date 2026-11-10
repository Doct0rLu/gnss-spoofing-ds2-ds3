function [score,gamma,alarm,joint,readyJoint]=deng_qcomposite(Q,valid,cleanQ,w,p)
%DENG_QCOMPOSITE 四个Q tap的和再平方；Q为历元×4×星。
% cleanQ同列顺序，按星冻结四点和方差（含协方差）；至少四星同时越阈。
[n,m,k]=size(Q);assert(m==4&&size(cleanQ,2)==4&&size(cleanQ,3)==k&&isequal(size(valid),[n k]));
assert(w>=1&&w==fix(w)&&p>0&&p<1);
u=reshape(sum(Q,2),n,k);s2=reshape(var(sum(cleanQ,2),1,1),1,k);
v=logical(valid)&isfinite(u)&isfinite(s2)&s2>0;
x=u.^2./s2;x(~v)=0;
score=movsum(x,[w-1 0],1);ready=movsum(double(v),[w-1 0],1)==w;score(~ready)=NaN;
gamma=2*gammaincinv(p,w/2,'upper');alarm=ready&score>gamma;
readyJoint=sum(ready,2)>=4;joint=readyJoint&sum(alarm,2)>=4;
end
