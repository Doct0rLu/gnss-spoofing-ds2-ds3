function path=prepare_tracking(project,receiver)
% 生成 CU-SDR 跟踪适配副本，保留原 GPL 归属，不修改第三方原件。
% 只移除图形进度条、增加批处理进度和积分结束复采样索引；不重写 DLL/PLL 环路。
% 每个替换目标要求唯一匹配，上游代码结构变化时立即失败，避免静默生成错误适配器。
% endSample 与 absoluteSample 配对，用于因果特征时间、连续性及离线相关器重算。
% 该函数处理的是 MATLAB 源码字符串；字符串内的英文属于上游内容和匹配锚点。
% Generate a local GPL-derived batch adapter; upstream remains untouched.
% Changes: headless progress; exact end-sample logging. DLL/PLL unchanged.
src=fileread(fullfile(receiver,'include','tracking.m'));
src=replace_once(src,'function [trackResults, channel]= tracking(fid, channel, settings)',...
    'function [trackResults, channel]= tracking_logged(fid, channel, settings)');
src=replace_once(src,'hwb = waitbar(0,''Tracking...'');','hwb = []; % Batch adapter: no GUI');
a=strfind(src,'CNoPos=get(hwb,'); b=strfind(src,'if (settings.fileType==1)');
assert(numel(a)==1 && numel(b)==1); src=[src(1:a-1) src(b:end)];
a=strfind(src,'            if (rem(loopCnt, 50) == 0)');
b=strfind(src,'            %% Read next block of data');
assert(numel(a)==1 && numel(b)==1);
progress=sprintf(['            if rem(loopCnt,5000)==0\n' ...
    '                fprintf(''Tracking PRN %%02d: %%d/%%d ms\\n'',channel(channelNr).PRN,loopCnt,codePeriods);\n' ...
    '            end\n\n']);
src=[src(1:a-1) progress src(b:end)];
src=replace_once(src,'trackResults.absoluteSample = zeros(1, settings.msToProcess);',...
    sprintf(['trackResults.absoluteSample = zeros(1, settings.msToProcess);\n' ...
    'trackResults.endSample = nan(1, settings.msToProcess);']));
src=replace_once(src,'            rawSignal = rawSignal'';',sprintf([...
    '            trackResults(channelNr).endSample(loopCnt) = ftell(fid)/dataAdaptCoeff/2;\n' ...
    '            rawSignal = rawSignal'';']));
src=replace_once(src,'                delete(hwb);','                % No GUI to delete.');
src=replace_once(src,'close(hwb)','% Batch adapter: no GUI to close.');
path=fullfile(project,'build','cusdr_gps_logged'); if ~isfolder(path), mkdir(path); end
fid=fopen(fullfile(path,'tracking_logged.m'),'wb'); assert(fid>0);
guard=onCleanup(@() fclose(fid)); fwrite(fid,src,'char');
end
function out=replace_once(src,a,b)
assert(numel(strfind(src,a))==1,'Unexpected upstream source: %s',a);
out=strrep(src,a,b);
end
