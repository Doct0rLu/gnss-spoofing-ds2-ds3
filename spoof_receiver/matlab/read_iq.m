function iq=read_iq(path,first,count)
% 读取 TEXBAT 原生复采样：小端 int16，交织顺序 I0,Q0,I1,Q1,…。
% first 是从零开始的复采样索引，count 是复采样数量；每个复采样占 4 字节。
% 返回 double 复数行向量；仅进行类型转换，不归一化、不重采样、不改变量化幅度。
% 显式检查边界和短读，避免把不完整文件或错误字节偏移当作合法信号。
% onCleanup 在正常返回和异常路径都关闭文件句柄。
% Zero-based complex-sample offset. TEXBAT native little-endian int16 IQ.
arguments
    path
    first (1,1) double {mustBeNonnegative,mustBeInteger}
    count (1,1) double {mustBeNonnegative,mustBeInteger}
end
assert(first+count<=flintmax/4,'Offset exceeds exact byte indexing');
info=dir(path); assert(~isempty(info) && mod(info.bytes,4)==0,'Invalid IQ file');
assert(4*(first+count)<=info.bytes,'IQ request beyond available samples');
fid=fopen(path,'rb','ieee-le'); assert(fid>0); guard=onCleanup(@() fclose(fid));
assert(fseek(fid,4*first,'bof')==0,'IQ seek failed');
[v,n]=fread(fid,[2,count],'int16=>double');
assert(n==2*count,'Short IQ read');
iq=complex(v(1,:),v(2,:));
end
