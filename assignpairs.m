function sites=assignpairs(sessnums,xinfo,plotparam)    
%assign single lfp pair to each recorded da site in xinfo compiled variable
subj='patra';
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites,subj);
uniquesites=unique({sites(1:end).site});
[dapair,lfppair]=getsitepairs(targdasites,subj);
for ises=1:length(sessnums)
    sid=find([sites.sessnum]==sessnums(ises));
    for ida=1:length(sid)
        daids=find(strcmp(dapair,sites(sid(ida)).probeid));
        targrow=find(contains({xinfo.siteda},sites(sid(ida)).probeid) & ... 
                strcmp({xinfo.sessionid},num2str(sessnums(ises))) &...
                contains({xinfo.sessiontype},'big')==1);
        curlfps=unique({xinfo(targrow).sitelfp});
        targlfps=lfppair(daids);
        flagempty=0;
        for iilfp=1:length(targlfps)
            lfpexist=find(strcmp(curlfps,targlfps{iilfp}));
            if ~isempty(lfpexist) && flagempty==0
                %restructure/replace da names in sites variable with lfp pair but
                %keep as separate identities
                lfpsitename=getlfpsites(sessnums(ises),targlfps(iilfp),subj);
                if ~isempty(lfpsitename) && ~(sid(ida)>1 && any(strcmp({sites(sid).probeid},lfpsitename.probeid)))
                    %if not bad site that day (e.g. 127, cl5)
                    %and not redundant for given session, ie, pair
                    %already assigend to another da site
                    sites(sid(ida)).daid=sites(sid(ida)).probeid;
                    sites(sid(ida)).dasite=sites(sid(ida)).site;
                    sites(sid(ida)).dach=sites(sid(ida)).ch;
                    sites(sid(ida)).probeid=lfpsitename.probeid;
                    sites(sid(ida)).site=lfpsitename.site;
                    sites(sid(ida)).ch=lfpsitename.ch;
                    flagempty=1;
                end
            end
        end 
        if isempty(sites(sid(ida)).dach)
            %still no site defined, just repeat then even though
            %redundant,     
            oids=1:length(sid);
            oids(ida)=[];
            cids=find(contains({sites(sid(oids)).probeid},'p'));
            if contains(sites(sid(ida)).probeid,'c')
            cids=find(contains({sites(sid(oids)).probeid},'c'));
            end
            lfpsite=sites(sid(cids(oids(1)))).site;
            lfpch=  sites(sid(cids(oids(1)))).ch;  
            lfpprobe=  sites(sid(cids(oids(1)))).probeid;  
            sites(sid(ida)).daid=sites(sid(ida)).probeid;
            sites(sid(ida)).dasite=sites(sid(ida)).site;
            sites(sid(ida)).dach=sites(sid(ida)).ch;
            sites(sid(ida)).probeid=lfpprobe;
            sites(sid(ida)).site=lfpsite;
            sites(sid(ida)).ch=lfpch;
            flagempty=1;
        end
    end        
end

end