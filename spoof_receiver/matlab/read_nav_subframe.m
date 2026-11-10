function frame=read_nav_subframe(prompt,start)
% Code-period index start may precede the record's first TLM. Missing bits
% remain '?' and are never supplied to bin2dec. HOW and payload must exist.
frame=struct('start_epoch',start,'bits',repmat('?',1,300),'available',false(1,10),...
    'parity_pass',false(1,10),'id',NaN,'tow_s',NaN,'tlm_verified',false);
for w=1:10
    first=start+600*(w-1)-40;last=first+639;
    if first<1 || last>numel(prompt),continue;end
    signs=sign(sum(reshape(prompt(first:last),20,32),1));signs(signs==0)=-1;
    raw=char('0'+(signs>0));
    frame.available(w)=true;frame.parity_pass(w)=navPartyChk(signs)~=0;
    if frame.parity_pass(w),frame.bits((w-1)*30+(1:30))=checkPhase(raw(3:32),raw(2));end
end
if frame.parity_pass(2)
    frame.id=bin2dec(frame.bits(50:52));frame.tow_s=mod(bin2dec(frame.bits(31:47))*6-6,604800);
end
frame.tlm_verified=frame.parity_pass(1) && strcmp(frame.bits(1:8),'10001011');
end
