function paths=prepare_texbat_iq(rawRoot,dataRoot)
%PREPARE_TEXBAT_IQ 从本地官方完整录制切出报告用的三个片段。
% rawRoot 包含 cleanStatic.bin、ds2.bin、ds3.bin，均从原录制 0 s 开始。
% dataRoot 省略时为本包 spoof_receiver/data。已有匹配片段会复用。
% 按原字节复制，不改变采样值或量化格式。文件匹配片段 SHA-256 后才改为正式名。
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'spoof_receiver','matlab'));
if nargin<2||isempty(dataRoot),dataRoot=fullfile(root,'spoof_receiver','data');end
if ~isfolder(dataRoot),mkdir(dataRoot);end
C=ds23_iq_contract();paths=cell(1,numel(C));
for k=1:numel(C)
 M=C(k);from=fullfile(rawRoot,M.official_file);to=fullfile(dataRoot,M.file);
 paths{k}=to;
 if isfile(to)
  f=dir(to);assert(f.bytes==M.expected_bytes&&strcmp(file_sha256(to),M.sha256), ...
   'Existing fragment does not match the report data: %s',to);
 else
  slice_iq_bytes(from,to,4*M.source_start_sample,M.expected_bytes,M.sha256);
 end
end
end
