function diagnose_ds23_prn7(tracks,s,M,t,score,g,out)
%DIAGNOSE_DS23_PRN7 重建弱峰图，并在不使用跟踪点的条件下搜索全码周期。
prns=[7 13 19 30]; C=cell(1,4);
for j=1:4
 tr=tracks([tracks.PRN]==prns(j));assert(isscalar(tr));
 C{j}=export_features(tr,s,M.source_start_sample);
end
f=C{1}; tr=tracks([tracks.PRN]==7);
amp=abs(complex(tr.I_P(:),tr.Q_P(:)))./(double(tr.endSample(:))-double(tr.absoluteSample(:)));
fig=figure('Visible','off','Color','w','Position',[50 50 970 845]);
guard=onCleanup(@() close(fig)); tl=tiledlayout(fig,4,1,'Padding','compact','TileSpacing','compact');
a=nexttile(tl);hold(a,'on');h=gobjects(4,1);color=lines(4);
for j=1:4
 F=C{j}; h(j)=plot(a,F.time_s(1:40:end),F.cn0_dbhz(1:40:end),'Color',color(j,:),'LineWidth',1.1);
end
ylabel(a,'C/N_0 (dB-Hz)');ylim(a,[10 66]);legend(a,h,compose('PRN %d',prns),'Location','northwest','NumColumns',4);
b=nexttile(tl);j=1:20:height(f);plot(b,f.time_s(j),f.phase_quality(j),'Color',color(1,:));
ylabel(b,'Phase quality');ylim(b,[-1 1.1]);
d=nexttile(tl);x=movmedian(20*log10(amp),200,'omitnan');plot(d,f.time_s(j),x(j));
ylabel(d,'Prompt amplitude (dB re 1 count)');
e=nexttile(tl);hold(e,'on');l=plot(e,t,score,'Color',[.75 .14 .11],'LineWidth',1.1);
y=yline(e,g,'k--');ylabel(e,'5 s Q-SQM');xlabel(e,'Original DS3 recording time (s)');
legend(e,[l y],{'Q-SQM','Clean 99th percentile'},'Location','northwest');
for ax=[a b d e],xlim(ax,[105 165]);grid(ax,'on');end
title(tl,'DS3 PRN 7 | Weak correlation peak and Q-SQM');
exportgraphics(fig,fullfile(out,'figures','diagnostic_prn07.png'),'Resolution',160);

fs=s.samplingFreq; N=25000; ca=generateCAcode(7); code=ca(ceil((1:N)*1.023e6/fs));
stamps=[110 115 119 121 125 130 135 140 145 149 155 160].';
freq=1000:100:2700; peak=nan(size(stamps)); doppler=peak; lag=peak;
for k=1:numel(stamps)
 raw=read_iq(s.fileName,stamps(k)*fs-M.source_start_sample-20*N,20*N);
 P=acquisition_power_map(raw,code,fs,freq);
 [v,i]=max(P(:)); [j,h]=ind2sub(size(P),i);
 doppler(k)=freq(j);lag(k)=h-1;peak(k)=10*log10(v/median(P(:)));
end
T=table(stamps,peak,doppler,lag,'VariableNames',{'window_end_s','peak_above_median_db','doppler_hz','lag_samples'});
writetable(T,fullfile(out,'prn7_independent_capture.csv'));
end
