function valid=nav_issue_consistent(clockPart,orbit2,orbit3)
% IS-GPS-200N: IODEs must equal each other and the low eight IODC bits.
values={clockPart.health,clockPart.IODC,orbit2.IODE_sf2,orbit3.IODE_sf3};
valid=all(cellfun(@(x) isscalar(x)&&isfinite(x),values));
if ~valid,return;end
valid=clockPart.health==0 && clockPart.IODC>=0 && clockPart.IODC<=1023 && ...
    mod(clockPart.IODC,256)==orbit2.IODE_sf2 && orbit2.IODE_sf2==orbit3.IODE_sf3;
end
