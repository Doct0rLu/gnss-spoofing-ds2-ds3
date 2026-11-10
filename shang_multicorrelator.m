function [score,gamma,alarm]=shang_multicorrelator(U,M,sigma2,p)
%SHANG_MULTICORRELATOR 复相关残差的白化能量；U每列为一个历元。
% 使用每实分量噪声方差sigma2；Prompt投影后共2(n-1)个自由度。
n=4*M+1;assert(size(U,1)==n&&M>=2&&M==fix(M)&&sigma2>0&&p>0&&p<1);
offset=(-2*M:2*M)'/M;h=max(1-abs(offset),0);
R=toeplitz(max(1-(0:n-1)/M,0));L=chol(R,'lower');
v=all(isfinite(U),1);score=nan(1,size(U,2));
z=L\(U(:,v)-h*U(2*M+1,v));score(v)=sum(abs(z).^2,1)/sigma2;
gamma=2*gammaincinv(p,n-1,'upper');alarm=isfinite(score)&score>gamma;
end
