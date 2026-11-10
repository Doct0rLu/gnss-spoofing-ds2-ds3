function run_ds23(stage,dataRoot,workRoot,workers,cleanReferenceFile)
%RUN_DS23 从官方 IQ 重建 DS2/DS3 报告数据。stage: all, receiver, derive。
% 从源码包根目录运行 run_ds23('all',dataRoot,workRoot)。
% dataRoot 需要三个片段：cleanStatic_first42s.bin、ds2_s85_245.bin、
% ds3_s85_245.bin。只有官方完整录制时，可先运行 prepare_texbat_iq。
% workers 默认为 0（串行）；大于 0 时需要 Parallel Computing Toolbox。
% 可选 cleanReferenceFile 为从 0 s 开始的官方 cleanStatic 完整录制，
% 仅用于重新生成图上配对干净峰标注，统计量和门限不使用它。
% 生成的图和绘图数据位于 workRoot/report。
% 示例：run_ds23('receiver',dataRoot,workRoot); run_ds23('derive',dataRoot,workRoot);
% 许可：help license_gpl_v2。原始数据来源：https://radionavlab.ae.utexas.edu/texbat/
root=fileparts(mfilename('fullpath')); project=fullfile(root,'spoof_receiver');
if nargin<1,stage='all';end
if nargin<2||isempty(dataRoot),dataRoot=fullfile(project,'data');end
if nargin<3||isempty(workRoot),workRoot=fullfile(project,'outputs');end
if nargin<4,workers=0;end
if nargin<5,cleanReferenceFile='';end
assert(ismember(string(stage),["all","receiver","derive"]),'Unknown stage');
addpath(root,fullfile(project,'matlab')); 
rx=fullfile(root,'third_party','CU-SDR-Collection','GPS','GPS_L1CA');
addpath(rx,fullfile(rx,'Common'),fullfile(rx,'include'),'-begin');
[ok,d]=fileattrib(dataRoot);assert(ok&&isfolder(dataRoot),'Supply the official IQ directory');dataRoot=d.Name;
if ~isfolder(workRoot),mkdir(workRoot);end
[~,w]=fileattrib(workRoot);workRoot=w.Name;
contract=ds23_iq_contract();
for k=1:3
 meta=contract(k);file=fullfile(dataRoot,meta.file);
 st=dir(file);assert(isscalar(st)&&st.bytes==meta.expected_bytes,'Missing or wrong size: %s',file);
 dest=fullfile(workRoot,meta.receiver_folder);
 if ~strcmp(stage,'derive')
  assert(strcmp(file_sha256(file),meta.sha256),'IQ SHA-256 mismatch: %s',file);
  run_ds23_receiver(file,meta,dest,project,rx,k==1);
 else
  J=load(fullfile(dest,'report.mat'),'report');J=J.report;
  assert(strcmp(J.source_iq_sha256,meta.sha256)&&J.source_start_sample==meta.source_start_sample, ...
      'Receiver outputs do not match the bundled IQ contract.');
  assert(strcmp(file_sha256(file),meta.sha256),'IQ SHA-256 mismatch: %s',file);
 end
end
if ~strcmp(stage,'receiver')
 build_ds23_data(dataRoot,workRoot,workers,cleanReferenceFile);
end
end
