function slice_iq_bytes(from,to,first,count,sha)
%SLICE_IQ_BYTES 逐块复制官方录制的指定字节区间，哈希吻合后再提交。
% first 是零基字节偏移，count 为字节数；不覆盖已有正式或临时片段。
assert(~isfile(to),'Destination already exists: %s',to);
part=[char(to) '.partial'];assert(~isfile(part),'Temporary fragment already exists: %s',part);
f=dir(from);assert(isscalar(f)&&f.bytes>=first+count,'Recording is missing or too short: %s',from);
a=fopen(from,'rb');assert(a>0,'Cannot open %s',from);in=onCleanup(@() fclose(a));
assert(fseek(a,first,'bof')==0,'Cannot seek to the requested source offset');
b=fopen(part,'wb');assert(b>0,'Cannot create %s',part);out=onCleanup(@() fclose(b));
remaining=count;
while remaining>0
 n=min(remaining,8*1024*1024);
 [bytes,got]=fread(a,n,'*uint8');assert(got==n,'Short read in source recording');
 assert(fwrite(b,bytes,'uint8')==n,'Short write in output fragment');
 remaining=remaining-n;
end
clear out in
assert(strcmp(file_sha256(part),sha),'Fragment SHA-256 did not match the report contract: %s',part);
[ok,message]=movefile(part,to);assert(ok,'%s',message);
end
