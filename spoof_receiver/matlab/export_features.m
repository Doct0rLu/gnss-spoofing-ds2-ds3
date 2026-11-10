function features=export_features(tracks,settings,sourceStart)
% 将跟踪状态转成检测器使用的特征表，不引入场景标签或外部真值。
% sourceStart 为片段在原始录制中的零基复采样偏移，用于统一到原录制时间轴。
% 每行对应完成的码周期，频率取本次相关实际使用的 NCO 值，不能替换为更新后的值。
% C/No 和相位质量只从已完成的历史窗口获得，再向后保持到下次更新，防止未来泄漏。
% completed 表示该积分数值完整；usable 额外要求预热、C/No 与相位质量通过工程门限。
% usable 只表示适合输出检测特征，不是信号真实性或锁定真实卫星的证明。
% 所有未完成/不可用历元仍导出；time_s 与 available_at_s 均使用原始录制积分结束时刻。
% One row per completed code period. Only past C/No/quality is available.
% NCO frequencies refer to the interval just correlated, before its update.
parts=cell(1,numel(tracks));
for c=1:numel(tracks)
    t=tracks(c); n=numel(t.I_P); ix=(1:n)';
    first=t.absoluteSample(:); last=t.endSample(:);
    ip=t.I_P(:); qp=t.Q_P(:); cn=nan(n,1); phaseQuality=nan(n,1);
    for k=1:numel(t.CNo.VSMIndex)
        at=t.CNo.VSMIndex(k);
        if at<40 || at>n, continue; end
        until=min(at+settings.CNo.VSMinterval-1,n);
        cn(at:until)=t.CNo.VSMValue(k);
        % 相位质量采用截至当前更新点的 40 个 Prompt 历元，不读取下一窗口。
        j=at-39:at; den=sum(ip(j).^2+qp(j).^2);
        if den>0, phaseQuality(at:until)=sum(ip(j).^2-qp(j).^2)/den; end
    end
    completed=isfinite(last) & last>first & isfinite(t.codeFreq(:)) & isfinite(t.carrFreq(:));
    % Explicit engineering availability gate, not proof of authentic lock.
    % 保持冻结质量门：400 历元预热、C/No≥20 dB-Hz、相位质量≥0.5。
    usable=completed & ix>400 & isfinite(cn) & cn>=20 & phaseQuality>=0.5;
    end_s=(last+sourceStart)/settings.samplingFreq;
    parts{c}=table(repmat(t.PRN,n,1),ix,uint64(first+sourceStart),uint64(last+sourceStart),...
        end_s,end_s,(last-first)/settings.samplingFreq,...
        t.I_E(:),ip,t.I_L(:),t.Q_E(:),qp,t.Q_L(:),...
        t.codeFreq(:),t.carrFreq(:)-settings.IF,cn,phaseQuality,completed,usable,...
        'VariableNames',{'prn','epoch','start_sample','end_sample','time_s','available_at_s',...
        'integration_s','i_e','i_p','i_l','q_e','q_p','q_l','code_hz','carrier_doppler_hz',...
        'cn0_dbhz','phase_quality','completed','usable'});
end
features=vertcat(parts{:});
end
