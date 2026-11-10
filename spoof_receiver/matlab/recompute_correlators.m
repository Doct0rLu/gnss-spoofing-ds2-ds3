function z=recompute_correlators(iqFile,track,epoch,settings,offsets,raw,code)
% 使用历史跟踪状态，在同一真实 IQ 历元上重算任意码偏移相关器。
% epoch 为跟踪历元索引，offsets 单位 chip；输出 z 为复相关累加值，顺序与 offsets 一致。
% 读取该历元实际样本范围，使用当时载波频率/残余相位及码 NCO，不运行或更新环路。
% 保留原幅度、数据位和相位，不做 Prompt 幅度归一化；归一化由明确的上层合同决定。
% 码相位采用与 CU-SDR 相同的冒号网格和 ceil 约定；浮点运算重排可能改变码边界采样。
% 批处理可传入本历元的原始 IQ 行向量 raw 和预生成的本星 C/A 码 code；省略则按原方式读取。
% 此接口支持相关形状展示和多相关器统计量，不能把跟踪 Prompt 当作独立真实峰。
% Offline diagnostic taps from the exact NCO states used for this epoch.
% Does not update loops, wipe bits or normalize amplitude.
first=track.absoluteSample(epoch); count=track.endSample(epoch)-first;
assert(isfinite(count) && count>0 && count==fix(count));
if nargin<6,raw=read_iq(iqFile,first,count);else,assert(isrow(raw)&&numel(raw)==count);end
if nargin<7,code=generateCAcode(track.PRN);end
t=(0:count-1)/settings.samplingFreq;
base=raw.*exp(-1i*(2*pi*track.carrFreq(epoch)*t+track.remCarrPhase(epoch)));
step=track.codeFreq(epoch)/settings.samplingFreq;
z=zeros(1,numel(offsets));
for j=1:numel(offsets)
    % Equivalent to CU-SDR's ceil(tcode)+1 on [last, code, first],
    % generalized periodically to arbitrary diagnostic code offsets.
    % Match upstream colon construction as well as its ceil convention;
    % differently associated floating-point operations near chip boundaries
    % can otherwise choose a neighboring code chip for one ADC sample.
    phase=(track.remCodePhase(epoch)+offsets(j)):step:...
        ((count-1)*step+track.remCodePhase(epoch)+offsets(j));
    assert(numel(phase)==count,'Code phase grid length mismatch');
    ix=mod(ceil(phase)-1,settings.codeLength)+1;
    z(j)=sum(code(ix).*base);
end
end
