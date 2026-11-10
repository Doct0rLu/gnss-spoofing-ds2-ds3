function z = recompute_correlators_vectorized(iqFile, track, epoch, settings, offsets, raw)
%RECOMPUTE_CORRELATORS_VECTORIZED 用一次码相位网格重算多个相关器。
% 输入 IQ 文件、CU-SDR 原始跟踪状态、历元索引、接收机设置和 chip 偏移；
% 输出 z 为与 offsets 同序的复相关累加值，不做 Prompt 或噪声归一化。
% 保持原始采样区间、载波 NCO、码 NCO、ceil 索引及逐 tap 求和顺序。
% 只复用同一历元的码相位网格；浮点边界若与参考实现不一致，必须
% 先由独立的逐样本对照核验，不能直接替换正式实验的参考路径。
first = track.absoluteSample(epoch);
count = track.endSample(epoch) - first;
assert(isfinite(count) && count > 0 && count == fix(count), ...
    '无效历元样本范围');
if nargin<6,raw=read_iq(iqFile,first,count);end
assert(isrow(raw)&&numel(raw)==count);
code = generateCAcode(track.PRN);
t = (0:count-1) / settings.samplingFreq;
base = raw .* exp(-1i * (2*pi*track.carrFreq(epoch)*t + ...
    track.remCarrPhase(epoch)));
step = track.codeFreq(epoch) / settings.samplingFreq;
startPhase = track.remCodePhase(epoch);
phase = startPhase:step:((count-1)*step + startPhase);
assert(numel(phase) == count, '码相位采样网格长度不符');
z = zeros(1, numel(offsets));
for j = 1:numel(offsets)
    index = mod(ceil(phase + offsets(j))-1, settings.codeLength) + 1;
    z(j) = sum(code(index) .* base);
end
end
