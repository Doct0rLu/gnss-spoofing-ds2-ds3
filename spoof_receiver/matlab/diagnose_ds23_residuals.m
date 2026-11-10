function diagnose_ds23_residuals(tracks,settings,meta,out)
%DIAGNOSE_DS23_RESIDUALS 分解典型时段的峰形残差与四点 Q 能量。
cases=[2 13 105;2 13 178;2 13 228;3 7 105;3 7 125;3 7 135;3 7 230;3 19 113;3 19 123;3 19 211];
values=nan(size(cases,1),7);
for j=1:size(cases,1)
 rec=cases(j,1);prn=cases(j,2);s=settings{rec};M=meta{rec};
 T=tracks{rec};tr=T([T.PRN]==prn);assert(isscalar(tr));
 t=(double(tr.endSample(:))+M.source_start_sample)/s.samplingFreq;
 idx=find(t<=cases(j,3),20,'last');assert(numel(idx)==20);
 U=zeros(81,numel(idx));
 for k=1:numel(idx)
  i=idx(k);n=tr.endSample(i)-tr.absoluteSample(i);
  U(:,k)=recompute_correlators_folded(s.fileName,tr,i,s).'/n;
 end
 ei=shang_multicorrelator(real(U),20,1,1e-4);
 eq=shang_multicorrelator(imag(U),20,1,1e-4);
 e=shang_multicorrelator(U,20,1,1e-4);
 u=sum(imag(U([31 36 46 51],:)),1);
 values(j,:)=[mean(e),mean(ei),mean(eq),mean(u.^2),mean(u)^2,var(u,1),mean(abs(U(41,:)).^2)];
end
T=array2table([cases values],'VariableNames',{'record_index','prn','end_s','residual','i_residual', ...
 'q_residual','q_sum_energy','q_sum_mean_squared','q_sum_variance','prompt_energy'});
writetable(T,fullfile(out,'correlator_cause.csv'));
end
