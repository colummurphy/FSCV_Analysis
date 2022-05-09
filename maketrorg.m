function trorg=maketrorg()
%make table organized for trial organization

trorg(1).history={'break','break','break'};
trorg(2).history={'small','small','small'};
trorg(3).history={'big','big','big'};
trorg(4).history={'break','break'};
trorg(5).history={'small','small'};
trorg(6).history={'big','big'};
trorg(7).history={'break'};
trorg(8).history={'small'};
trorg(9).history={'big'};
trorg(10).future={'reward'};
trorg(11).future={'break'};
trorg(12).future={'reward','reward'};
trorg(13).future={'break','break'};
trorg(14).future={'reward','reward','reward'};
trorg(15).future={'break','break','break'};
trorg(16).history={'break','break','break','break'};
trorg(17).history={'small','small','small','small'};
trorg(18).history={'big','big','big','big'};

for ii=1:length(trorg)
    hlab='';
    if ~isempty(trorg(ii).history)
        hlab=trorg(ii).history{1};
        nums=find(contains(trorg(ii).history,hlab));
        trorg(ii).label=['-' num2str(length(nums)) hlab];
    end
    flab='';
    if ~isempty(trorg(ii).future)
        flab=trorg(ii).future{1};
        nums=find(contains(trorg(ii).future,flab));
        trorg(ii).label=['+' num2str(length(nums)) flab];
    end
end

end