function [eph,subFrameStart,TOW,audit]=NAVdecoding_checked(prompt,settings)
% Offline LNAV validation. No PRN/scenario/location/ephemeris input.
% Search rule follows CU-SDR; word parity covers all consumed fields.
eph=eph_structure_init();subFrameStart=Inf;TOW=Inf;
audit=struct('status','no_valid_anchor','anchor_epoch',NaN,'anchor_tow_s',NaN,...
    'selected_start_epochs',nan(1,3),'selected_partial_tlm',false(1,3),...
    'selected_issue',NaN,'frames',{{}},'available_after_epoch',NaN);
bits=ones(size(prompt));bits(prompt<=0)=-1;
preamble=kron([1 -1 -1 -1 1 -1 1 1],ones(1,20));cor=xcorr(bits,preamble);
middle=(numel(cor)+1)/2;
candidates=find(abs(cor(middle:2*middle-1))>153);
candidates=candidates(candidates>40 & candidates+1199<=numel(prompt));
for index=reshape(candidates,1,[])
    if ~any(candidates==index+6000),continue;end
    frame=read_nav_subframe(prompt,index);
    if frame.tlm_verified && frame.parity_pass(2) && ismember(frame.id,1:5)
        subFrameStart=index;TOW=frame.tow_s;anchorID=frame.id;break;
    end
end
if ~isfinite(subFrameStart),return;end
audit.anchor_epoch=subFrameStart;audit.anchor_tow_s=TOW;
parts=cell(1,3);selected=[];frames={};
firstShift=floor((1-subFrameStart)/6000)-1;lastShift=ceil((numel(prompt)-subFrameStart)/6000);
for shift=firstShift:lastShift
    start=subFrameStart+6000*shift;
    if start+560<1 || start+5999>numel(prompt),continue;end
    frame=read_nav_subframe(prompt,start);
    expectedID=mod(anchorID-1+shift,5)+1;expectedTOW=mod(TOW+6*shift,604800);
    consistent=frame.parity_pass(2) && frame.id==expectedID && frame.tow_s==expectedTOW;
    partialTLM=~frame.available(1) && start<41;
    payloadOK=consistent && (frame.tlm_verified || partialTLM) && all(frame.parity_pass(3:10));
    reason='accepted_payload';
    if ~consistent,reason='how_or_sequence_invalid';
    elseif ~(frame.tlm_verified || partialTLM),reason='tlm_invalid';
    elseif ~all(frame.parity_pass(3:10)),reason='payload_parity_failed';end
    item=struct('start_epoch',start,'subframe_id',frame.id,'tow_s',frame.tow_s,...
        'available_words',find(frame.available),'failed_words',find(frame.available & ~frame.parity_pass),...
        'partial_tlm',partialTLM,'payload_accepted',payloadOK,'reason',reason);
    frames{end+1}=item; %#ok<AGROW>
    if ~payloadOK || ~ismember(frame.id,1:3),continue;end
    part=struct('eph',decode_nav_subframe_checked(eph_structure_init(),frame.bits),...
        'start',start,'partial_tlm',partialTLM,'bits',frame.bits);
    parts{frame.id}=part;
    if isempty(selected) && all(~cellfun(@isempty,parts))
        a=parts{1}.eph;b=parts{2}.eph;c=parts{3}.eph;
        if nav_issue_consistent(a,b,c)
            selected=parts;
        end
    end
end
audit.frames=frames;
if isempty(selected)
    audit.status='no_consistent_checked_ephemeris';eph=eph_structure_init();return;
end
for k=1:3
    eph=decode_nav_subframe_checked(eph,selected{k}.bits);
    audit.selected_start_epochs(k)=selected{k}.start;
    audit.selected_partial_tlm(k)=selected{k}.partial_tlm;
end
eph.TOW=TOW;audit.status='validated';audit.selected_issue=mod(eph.IODC,256);
audit.available_after_epoch=max(audit.selected_start_epochs)+5999;
end
