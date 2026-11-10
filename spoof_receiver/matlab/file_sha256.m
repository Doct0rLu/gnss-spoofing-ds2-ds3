function hash=file_sha256(path)
% Stream a large IQ file; memory does not depend on total recording length.
md=java.security.MessageDigest.getInstance('SHA-256');
fid=fopen(path,'rb'); assert(fid>0); guard=onCleanup(@() fclose(fid));
while true
    v=fread(fid,8e6,'*uint8'); if isempty(v), break; end
    md.update(typecast(v,'int8'));
end
hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2)',1,[]));
end
