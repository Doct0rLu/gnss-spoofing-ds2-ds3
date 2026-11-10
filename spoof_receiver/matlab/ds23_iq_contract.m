function C=ds23_iq_contract()
%DS23_IQ_CONTRACT 报告使用的官方 TEXBAT 片段、采样起点及 SHA-256。
% 25 MS/s，零中频，小端有符号 int16，I/Q 交织，每复采样 4 字节。
whole={'cleanStatic.bin','ds2.bin','ds3.bin'};
files={'cleanStatic_first42s.bin','ds2_s85_245.bin','ds3_s85_245.bin'};
rx={'clean_full','ds2_pull','ds3_pull'};
start=[0 85 85];duration=[42 160 160];fs=25e6;
sha={'34cc75398d4d9ff73f81bec77f182d585597a76f726679e48d73e1bd105b3a4b', ...
     'f337a60f4f87138cf14155fe7cbd4fc3699f7a83a4173b13036e4777c9059016', ...
     'ba18a9b7a401ccafbdec9a261cda08f9fd913b27ed75c052b6af4273e7131047'};
for k=1:numel(files)
 C(k)=struct('official_file',whole{k},'file',files{k},'receiver_folder',rx{k}, ...
  'source_start_sample',start(k)*fs,'sample_count',duration(k)*fs,'sample_rate_hz',fs, ...
  'if_hz',0,'expected_bytes',duration(k)*fs*4,'sha256',sha{k},'status','complete'); %#ok<AGROW>
end
end
