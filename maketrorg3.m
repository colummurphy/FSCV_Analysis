function trorg=maketrorg3()
%make table organized for trial organization
trorg={};
trorg(1).history={'target','target','target'};
trorg(2).history={'fix','fix','fix'};

trorg(3).history={'small','small','small'};
trorg(4).history={'big','big','big'};
trorg(5).history={'target','target'};
trorg(6).history={'fix','fix'};

trorg(7).history={'small','small'};
trorg(8).history={'big','big'};
trorg(9).history={'fix'};
trorg(10).history={'target'};

trorg(11).history={'small'};
trorg(12).history={'big'};
trorg(13).future={'reward'};
trorg(14).future={'target'};
trorg(15).future={'fix'};

trorg(16).future={'reward','reward'};
trorg(17).future={'fix','fix'};
trorg(18).future={'target','target'};
trorg(19).future={'reward','reward','reward'};
trorg(20).future={'fix','fix','fix'};
trorg(21).future={'target','target','target'};
trorg(22).history={'target','target','target','target'};
trorg(23).history={'fix','fix','fix','fix'};
trorg(24).history={'small','small','small','small'};
trorg(25).history={'big','big','big','big'};

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