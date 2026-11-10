function P=acquisition_power_map(raw,code,fs,freqs)
% 独立捕获功率面：以码延迟和搜索频率为二维坐标。
% raw 为复数 IQ，code 为一个 1 ms 的采样本地码，fs 为采样率 Hz，freqs 为搜索频率 Hz。
% 输出 P 的行对应 freqs，列对应循环码延迟采样点；绝对单位为相关幅度平方。
% 每个 1 ms 块先去载波，再做 FFT 循环相关，最后跨块非相干累加功率。
% 跨毫秒累加功率可避免导航位翻转造成整段相干抵消；不足 1 ms 的末尾不参与。
% 该图只呈现信号结构；峰值不是经过认证的真实/欺骗标签，也不输入检测门限。
% Circular 1ms FFT acquisition, summed NONcoherent power across whole ms.
% No NCO history, tracking Prompt, scenario labels or truth used here.
N=round(fs/1000);K=floor(numel(raw)/N);assert(K>=1 && numel(code)==N);
blocks=reshape(raw(1:N*K),N,K);time=(0:N-1)'/fs;replica=conj(fft(code(:)));
P=zeros(numel(freqs),N);
for j=1:numel(freqs)
    mixed=blocks.*exp(-1i*2*pi*freqs(j)*time);
    % FFT 乘以本地码频谱共轭后 IFFT，实现所有循环码延迟的并行相关。
    correlation=ifft(fft(mixed,[],1).*replica,[],1);
    P(j,:)=sum(abs(correlation).^2,2)';
end
end
