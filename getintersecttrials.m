function trialgrps=getintersecttrials(trialnums,intersectgrps,varargin)
%intersect trials across groups of select groups of trialtype-grouped trialnums
%trialnums{1} may include trialnums{grp1}(type1) &
%trialnums{grp1}(type2), trialnums{grp2}(type1) ...etc.
%intersect {1} grp1 / type1 with {2} grp1 / type 1 & {1} grp2 / type1 with {2} grp 2 / type 1 etc.
%output trialnums with separate groups and types as output by getmulttrials
%each group should have same # of grps and same # of types to intersect
%across
trialgrps=[];
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dapos'
                          
    end
    argnum=argnum+1;
end

for ix=1:length(intersectgrps)
    trialgrps(ix).site=[];
    trialgrps(ix).cat=[];
    for icat=1:length(intersectgrps{ix})
        xcat=intersectgrps{ix}{icat};
        targrow=find(contains({trialnums.cat},xcat));
        if icat==1
            trialgrps(ix).site=trialnums(targrow).site;
            trialgrps(ix).cat=trialnums(targrow).cat;
        else
            trialgrps(ix).site=[trialgrps(ix).site 'and' trialnums(targrow).site];
            trialgrps(ix).cat=[trialgrps(ix).cat 'and' trialnums(targrow).cat];
        end
        for itype=1:length(trialnums(targrow).trialgrps)
            if icat==1
                trialgrps(ix).trialgrps(itype).trials=trialnums(targrow).trialgrps(itype).trials;
                trialgrps(ix).trialgrps(itype).type=trialnums(targrow).trialgrps(itype).type;
            else
                trialgrps(ix).trialgrps(itype).trials=intersect(trialgrps(ix).trialgrps(itype).trials, trialnums(targrow).trialgrps(itype).trials);
            end
        end
    end
end
    
end