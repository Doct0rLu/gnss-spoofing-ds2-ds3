function [score,gamma,alarm,joint,readyJoint,accel]=feng_cd2p(time,rho,F,sigma2,L,p,noiseModel)
%FENG_CD2P 规则时间上的伪距二阶差分，扣除已知卫星项F后累积。
% rho为米，F为m/s^2；sigma2为冻结的逐星方差。F=0是显式的工程近似。
% 默认沿用能量和；iid_range 模式对独立原伪距噪声的重叠二阶差分白化。
if nargin<7,noiseModel='nominal';end
time=time(:);[n,k]=size(rho);dt=median(diff(time));
assert(n==numel(time)&&isequal(size(F),size(rho))&&numel(sigma2)==k&&dt>0&&L==fix(L)&&L>=1&&p>0&&p<1);
good=abs(diff(time)-dt)<=max(1e-6,1e-3*dt);accel=nan(n,k);
accel(3:end,:)=diff(rho,2,1)/dt^2-F(3:end,:);
accel(~[false;false;good(1:end-1)&good(2:end)],:)=NaN;
s2=reshape(sigma2,1,[]);v=isfinite(accel)&isfinite(s2)&s2>0;
x=accel.^2./s2;x(~v)=0;ready=movsum(double(v),[L-1 0],1)==L;
score=movsum(x,[L-1 0],1);score(~ready)=NaN;
if strcmp(noiseModel,'iid_range')
 c=zeros(1,L);r=[1 -2/3 1/6];c(1:min(L,3))=r(1:min(L,3));
 B=chol(toeplitz(c),'lower');
 for j=1:k
  for i=find(ready(:,j)).'
   a=B\accel(i-L+1:i,j);score(i,j)=sum(a.^2)/s2(j);
  end
 end
else,assert(strcmp(noiseModel,'nominal'));end
gamma=2*gammaincinv(p,L/2,'upper');alarm=ready&score>gamma;
visible=isfinite(rho);nv=sum(visible,2);readyJoint=nv>0&all(ready|~visible,2);
joint=readyJoint&sum(alarm&visible,2)>=nv/2;
end
