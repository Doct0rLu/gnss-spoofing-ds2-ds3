function verify_source_smoke(dataRoot,trackRoot,referenceData)
%VERIFY_SOURCE_SMOKE 从原始 DS3 IQ 复算少量历元，核对弱峰高台与冻结绘图数据。
% trackRoot 是用 run_ds23('receiver',...) 得到的目录。
% 可选 referenceData 是轻量重绘包的 figure_data.mat，用来逐点对照热图。
root=fileparts(mfilename('fullpath'));
project=fullfile(root,'spoof_receiver');
rx=fullfile(root,'third_party','CU-SDR-Collection','GPS','GPS_L1CA');
addpath(root,fullfile(project,'matlab'),rx,fullfile(rx,'Common'),fullfile(rx,'include'));
C=ds23_iq_contract();M=C(3);
R=load(fullfile(trackRoot,'ds3_pull','tracking.mat'),'tracks','settings');
s=R.settings;s.fileName=fullfile(dataRoot,'ds3_s85_245.bin');tr=R.tracks([R.tracks.PRN]==7);
A=derive_ds23_star(tr,s,M,[132 132.16],true,true);
assert(numel(A.time)>100&&all(A.valid)&&all(isfinite(A.energy)));
[peak,k]=max(vecnorm(A.sunQ,2,2)./abs(A.promptI));
assert(abs(A.time(k)-132.148)<.001&&abs(peak/5000-23.42)<.08, ...
    'The PRN 7 near-zero Prompt-I event did not match the report');
assert(isequal(A.capture.t,132)&&all(isfinite(A.capture.I),'all')&&all(isfinite(A.capture.Q),'all'));
fprintf('PRN 7 high plateau contribution %.5f; DS3 81-tap I/Q generated at 132 s.\n',peak/5000);
if nargin>=3&&~isempty(referenceData)
 S=load(referenceData,'D');
 D=S.D(strcmp({S.D.id},'sun_ds3_prn07'));i=find(D.capture.t==132);
 assert(isscalar(i)&&max(abs(double(A.capture.I)-double(D.capture.I(i,:))))<.001);
 assert(max(abs(double(A.capture.Q)-double(D.capture.Q(i,:))))<.001);
 fprintf('Both 81-tap branches match the frozen report window.\n');
end
end
