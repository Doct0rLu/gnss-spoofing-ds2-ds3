function batch=prepare_navigation(project,receiver)
% GPL-derived adapters; preserve upstream source and copyright headers.
batch=fullfile(project,'build','cusdr_nav_checked');if ~isfolder(batch),mkdir(batch);end
src=fileread(fullfile(receiver,'include','ephemeris.m'));
header=src(1:strfind(src,'%% Check if there is enough data')-1);
header=strrep(header,'function [eph, TOW] = ephemeris(bits, D30Star)',...
    'function eph=decode_nav_subframe_checked(eph,subframe)');
a=strfind(src,'    switch subframeID');b=strfind(src,'    end % switch subframeID');
assert(numel(a)==1 && numel(b)==1);
body=src(a:b+length('    end % switch subframeID')-1);
% IS-GPS-200N 20.3.3.3.1.5: IODC low bits are word8 bits1:8 (211:218).
body=strrep(body,'[subframe(83:84) subframe(197:204)]','[subframe(83:84) subframe(211:218)]');
body=strrep(body,'eph.idValid(1)  = 2;','eph.idValid(2)  = 1;');
body=strrep(body,'eph.idValid(3)  = 3;','eph.idValid(3)  = 1;');
prefix=sprintf(['\n%% Local parser: fields only from parity-checked words.\n' ...
    'assert(ischar(subframe) && numel(subframe)==300);\n' ...
    'assert(all(subframe(31:300)==''0'' | subframe(31:300)==''1''),''Missing required HOW/payload bits'');\n' ...
    'gpsPi=3.1415926535898; subframeID=bin2dec(subframe(50:52));\n']);
write_source(fullfile(batch,'decode_nav_subframe_checked.m'),[header prefix body newline 'end' newline]);
src=fileread(fullfile(receiver,'include','postNavigation.m'));
src=strrep(src,'function [navSolutions, eph] = postNavigation(trackResults, settings)',...
    'function [navSolutions, eph, navAudit] = postNavigation_checked(trackResults, settings)');
src=strrep(src,'if (settings.msToProcess < 36000)',sprintf('navAudit=cell(1,settings.numberOfChannels);\nif (settings.msToProcess < 36000)'));
src=strrep(src,'[eph(PRN), subFrameStart(channelNr), TOW(channelNr)]',...
    '[eph(PRN), subFrameStart(channelNr), TOW(channelNr), navAudit{channelNr}]');
src=strrep(src,'NAVdecoding(trackResults(channelNr).I_P,settings);',...
    'NAVdecoding_checked(trackResults(channelNr).I_P,settings); navAudit{channelNr}.prn=PRN;');
write_source(fullfile(batch,'postNavigation_checked.m'),src);
end
function write_source(path,src)
fid=fopen(path,'wb');assert(fid>0);guard=onCleanup(@() fclose(fid));fwrite(fid,src,'char');
end
