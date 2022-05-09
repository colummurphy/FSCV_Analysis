function trorg=maketrorg4()
%make table organized for trial organization

trorg(1).history={'break','break','break'};
trorg(2).history={'reward','reward','reward'};
trorg(3).history={'break','break'};
trorg(4).history={'reward','reward'};
trorg(5).history={'break'};
trorg(6).history={'reward'};
trorg(9).future={'reward'};
trorg(10).future={'break'};
trorg(11).future={'reward','reward'};
trorg(12).future={'break','break'};
trorg(13).future={'reward','reward','reward'};
trorg(14).future={'break','break','break'};
trorg(7).history={'reward','reward','reward','reward'};
trorg(8).history={'break','break','break','break'};

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