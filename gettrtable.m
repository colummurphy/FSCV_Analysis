function trtable=gettrtable(trlists,trorg)
%get trial nums for current 'big','small',etc. based on previous trial
%history or trial future as supplied in trorg for al sessions in trlists

trtable={};
for ii=1:length(trlists)
    list=trlists(ii).list;
    types=unique({list.type});
    trtable(ii).sessid=trlists(ii).sessid;
    for it=1:length(types)
        trtable(ii).trorg(it).type=types{it};
        trseq={list.type};        %trial sequence types
        trseqtnums=[list.id];       %trial nums for each
        for oo=1:length(trorg)
            hname='';
            if ~isempty(trorg(oo).history)
                hids={};
                hc=0;
                hinall=[];      %overlap with other trials
                for ih=1:length(trorg(oo).history)
                    hname=trorg(oo).history{ih};
                    hids{ih}=find(contains(trseq,hname))-hc;        %subtract index to check if also previous trial name included, 
                    hc=hc+1;
                    if ih>1
                    hinall=intersect(hinall,hids{ih});
                    end
                    if ih==1
                        hinall=hids{ih};
                    end
                end
                %includ ecurrent trial type in intersection
                hids{length(trorg(oo).history)+1}=find(contains(trseq,types{it}))-hc; %subtract count to make sure current trial type after targeted trials
                hinall=intersect(hinall,hids{length(trorg(oo).history)+1});
                confirmcurtypes=all(contains(trseq(hinall+hc),types{it})); %make sure all identified types correspond to current trail type;
                if ~confirmcurtypes
                    warning(['sess ' trtable(ii).sessid ' type ' types{it} ' not all current types, need to check' ])
                end
                targetts=trseqtnums(hinall+hc)-99;        %default, 100 first trial
                trtable(ii).trorg(it).grp(oo).trials=targetts;
                trtable(ii).trorg(it).grp(oo).label=trorg(oo).label;
            end
            fname='';
            if ~isempty(trorg(oo).future)
                hids={};
                hc=0;
                hinall=[];      %overlap with other trials
                for ih=1:length(trorg(oo).future)
                    fname=trorg(oo).future{ih};
                    hids{ih}=find(contains(trseq,fname))-hc;        %subtract index to check if also previous trial name included, 
                    hc=hc+1;
                    if ih>1
                    hinall=intersect(hinall,hids{ih});
                    end
                    if ih==1
                        hinall=hids{ih};
                    end
                end
                %includ ecurrent trial type in intersection
                hids{length(trorg(oo).future)+1}=find(contains(trseq,types{it}))+1; %here add only 1 to make sure before first identified trial in those showing seq
                hinall=intersect(hinall,hids{length(trorg(oo).future)+1});
                confirmcurtypes=all(contains(trseq(hinall-1),types{it})); %make sure all identified types correspond to current trail type;
                if ~confirmcurtypes
                    warning(['sess ' trtable(ii).sessid ' type ' types{it} ' not all current types, need to check' ])
                end
                targetts=trseqtnums(hinall-1)-99;        %default, 100 first trial
                trtable(ii).trorg(it).grp(oo).trials=targetts;
                trtable(ii).trorg(it).grp(oo).label=trorg(oo).label;
            end
        end
    end
end
