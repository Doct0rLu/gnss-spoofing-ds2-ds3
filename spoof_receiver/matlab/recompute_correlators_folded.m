function z=recompute_correlators_folded(iqFile,tr,k,settings,raw)
%RECOMPUTE_CORRELATORS_FOLDED 精确码相位分箱加速固定81点相关，不重采样IQ。
% 输入仍为原始IQ和历史NCO；输出-2:0.05:2 chip的复数累加值。
% 对整数tap j，ceil(phi+j/20)=ceil((ceil(20phi)+j)/20)。把载波
% 擦除后的采样按ceil(20phi)折叠到20460格，再作循环相关；这是求和
% 重排，不改变相关模板。浮点边界附近回退逐tap参考，以保持ceil约定。
% 必须通过verify_fullrate_v09的全部真实历元/判决核验后才用于新实验。
persistent cachedPrn cachedCodeFft
n=tr.endSample(k)-tr.absoluteSample(k);offsets=(-40:40)/20;
if nargin<5,raw=read_iq(iqFile,tr.absoluteSample(k),n);end
assert(isrow(raw)&&numel(raw)==n);
step=tr.codeFreq(k)/settings.samplingFreq;
phase=tr.remCodePhase(k):step:((n-1)*step+tr.remCodePhase(k));
assert(numel(phase)==n);
if any(abs(20*phase-round(20*phase))<1e-8)
    z=recompute_correlators_vectorized(iqFile,tr,k,settings,offsets,raw);return;
end
if isempty(cachedPrn)||cachedPrn~=tr.PRN
    code=generateCAcode(tr.PRN);grid=0:20459;
    template=code(mod(ceil(grid/20)-1,1023)+1);
    cachedCodeFft=conj(fft(template));cachedPrn=tr.PRN;
end
base=raw.*exp(-1i*(2*pi*tr.carrFreq(k)*(0:n-1)/settings.samplingFreq+tr.remCarrPhase(k)));
bins=accumarray((mod(ceil(20*phase),20460)+1).',base.',[20460 1],@sum,complex(0));
c=ifft(fft(bins.').*cachedCodeFft);
z=c(mod(-(-40:40),20460)+1);
end
