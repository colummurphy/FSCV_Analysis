function pct = getpctdirect(isub, varargin)
%01/29/2019 include M component, take out bg component
%Use patra_clean_3pcs_nobg_frommovep2_orig4182018.mat for Ats_mat and
%Cts_mat
%Originally was using patra_clean_titanium_da_ph_bg_movep2.mat
%NEED TO UPDATE CLEO TEMPALTES TO HAVE THIS STRUCTURE (IE NO BG)
%PREVIOUS SCRIPT SAVEDS AS getpctdirectold 01/28/2019
%get pca directly from isub, and already found glitches (input argument)
smoothwidth=4;        %# of samples +/- to interpolate with nan around glitches/artifacts
argnum=1;
qthres=[];
ats=[];
cts=[];
includem=1;         %01/28/2019
rthresm=0.85;
rthresp=0.95;    %rthresout for ph
rthresb=0.95;  %rthresout for bg
imaxbg=[];
glitchids=[];
numpcs=3;       %da,ph,M pc components to use default 01/2019
rthresb=0.85;  %rthresout for bg 1/2019
mstdlim=2; %1/2019

while argnum <= length(varargin)
    switch varargin{argnum}
        case 'qthres'
            %qthres (calculated by pcr of just dopamine/ph components)
            argnum=argnum+1;
            qthres=varargin{argnum};
        case 'ats'
            %CV templates for pca
            argnum=argnum+1;
            ats=varargin{argnum};            
        case 'cts'
            %conc scale values for pca
            argnum=argnum+1;
            cts=varargin{argnum};  
        case 'nom'
            %do not include M (K(:,3)) in calculating dopamine or not,
            %1/2019
            includem=0;
        case 'rthresm'
            %movement component correlation coefficient threshold exclusion
            argnum=argnum+1;
            rthresm=varargin{argnum};
        case 'rthresp'
            %ph component correlation coefficient threshold exclusion
            argnum=argnum+1;
            rthresp=varargin{argnum};
        case 'rthresb'
            %bg component correlation coefficient threshold exclusion
            argnum=argnum+1;
            rthresb=varargin{argnum};
        case 'imaxbg'
            %max absolute current for bg current for scaling concentration
            argnum=argnum+1;
            imaxbg=varargin{argnum};
        case 'glitchids'
            %glitch ids to nan (detected in findglitches script)
            argnum=argnum+1;
            glitchids=varargin{argnum};
        case 'nanwidth'
            %width of nanned signal padding
            argnum=argnum+1;
            smoothwidth=varargin{argnum};
        case 'numpcs'
            %num pcs to use, ie. da, ph, bg, m
            argnum=argnum+1;
            numpcs=varargin{argnum};

        otherwise
            error('badoption');
    end
    argnum = argnum + 1;
end
    
idsm=[];        %idx of signals  correlated to movement template ie. K(:,4)
KBG=[];
KM=[];
KPH=[];
KscaledM=[];        %movement templates repeated to find correlated data (ie. artifacts)
Mplot=[];           %also calculate M current component if we were to include all 4 signal categories
%if do not want to use movement component, but exists in loaded pca
%templates, remove from templates cts and ats variables
%still calculated KscaledM for artifact removal
%{
if includem==0 && size(cts,1)>=4
    [Vc, F, ~, ~]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %current contributions predicted of each component 
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    Mplot=[];
    KM=zeros(length(Vc),1);
    KscaledM=KM*Dproj(3,:);
    try
        Mplot=Cu(4,:);
        KM=K(:,4);
        KscaledM=KM*Dproj(4,:);     %ideal CV template to correlate to
    catch
    end
    KscaledM=KM*Dproj(4,:);     %ideal CV template to correlate to
    
    %now reduce templates provided to without movement templates from 
    %predict concentrations of main non-overlapping components
    midx=find(cts(4,:)==0);
    midx=midx(end);
    cts=cts(1:3,1:midx);
    ats=ats(:,1:midx);
elseif includem==1 
    [Vc, F, ~, ~]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %current contributions predicted of each component 
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    Mplot=[];
    KM=zeros(length(Vc),1);
    KscaledM=KM*Dproj(3,:);    
end
%}
dsdata=[];
K=[];

if numpcs==3
    %REMOVE LAST COMPONENT< BG, 1 /2019 after getting BG component for
%reference
if size(cts,1)==4
[Vc, F, Aproj, Qa]=generatePCR(ats,cts);
Dproj=Vc'*isub;
Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
E=isub-(Vc*Dproj);     %residuals
Q=diag(E'*E)';      %nA^2
BGplot=Cu(4,:);          %1/2019, changed from BGplot
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
KBG=K(:,4);
end
if size(cts,1)>3
    midx=find(cts(4,:)==0);
    midx=midx(end);
    cts=cts(1:3,1:midx);
    ats=ats(:,1:midx);
end
[Vc, F, Aproj, Qa]=generatePCR(ats,cts);
Dproj=Vc'*isub;
Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
E=isub-(Vc*Dproj);     %residuals
Q=diag(E'*E)';      %nA^2
IDAproj=Cu(1,:);            %current da projection
pHplot=Cu(2,:);
Mplot=Cu(3,:);          %1/2019, changed from BGplot
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
KDA=K(:,1);
KPH=K(:,2);
KM=K(:,3);          %3RD COMPONENT NOW MOVEMENT, NOT BG 01/2019
KscaledM=KM*Dproj(3,:);     %ideal CV template to correlate to

elseif numpcs==4
    %DA,PH,M,BG
        [Vc, F, Aproj, Qa]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
    E=isub-(Vc*Dproj);     %residuals
    Q=diag(E'*E)';      %nA^2
    IDAproj=Cu(1,:);            %current da projection
    pHplot=Cu(2,:);
    Mplot=Cu(3,:);
    BGplot=Cu(4,:);     %make bg zeros
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    KDA=K(:,1);
    KPH=K(:,2);
    KM=K(:,3);
    KBG=K(:,4);
    KscaledM=KM*Dproj(3,:);     %ideal CV template to correlate to
elseif numpcs==2
    %reduce to 2 components only
    if size(cts,1)>2
        midx=find(cts(3,:)~=0);
        midx=midx(1)-1;
        cts=cts(1:2,1:midx);
        ats=ats(:,1:midx);
    end
       
    [Vc, F, Aproj, Qa]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
    E=isub-(Vc*Dproj);     %residuals
    Q=diag(E'*E)';      %nA^2
    IDAproj=Cu(1,:);            %current da projection
    pHplot=Cu(2,:);
    BGplot=nan(1,size(Cu,2));     %make bg zeros
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    KDA=K(:,1);
    KPH=K(:,2);
    KBG=nan(size(K,1),1);
end

[aa,oxid]=max(K(:,1));   %DA oxidation current


%get correlation coefficients between data all time and movementtemplate
if exist('smoothdata')
    %if is function
    dsdata=smoothdata(isub,2,'gaussian',smoothwidth); 
else
    %use <2013 version of smooth data, need to check more carefully
    dsdata=smoothts(isub,'g',smoothwidth); 
end

rr=[];
rm=[];
rb=[];
routb=[];
rp=[];
routp=[];
if ~isempty(KscaledM)
rr=abs(corr(dsdata,KscaledM(:,1)));
rm=find(rr>rthresm);
end
if ~isempty(KBG)
rb=abs(corr(dsdata,KBG));
routb=find(rb>rthresb);
end

rp=abs(corr(dsdata,KPH));
routp=find(rp>rthresp);

if isempty(idsm)
    idsm=rm';
else
    if ~isempty(rm)
    idsm=[idsm rm'];        %idx of bad signals (ie. correlated to movement template or glitch)
    end
end

idsm=[idsm routb' routp'];

% ITplot=ITplot./(0.0341*ImaxBG)*1000;    %normalize to conc
ITplot=IDAproj./(0.0282*imaxbg*1.24)*1000;    %normalize to conc
DAproj=ITplot;

idsqout=find(Q>(qthres));     %Qa forDA & PH
idsout=[idsm glitchids idsqout];     %if above Rthres then must be artifact
ITplot(idsout)=nan;    
Mthres=nanstd([pHplot IDAproj BGplot])*mstdlim;     %1/20/2019 Thres for value of M component (not just corr looking for);
IDsaboveMthres=find(Mplot>Mthres);
ITplot(IDsaboveMthres)=nan;

pct.pHplot=pHplot;         %current's before nan
pct.Mplot=Mplot;          %current's before nan
pct.BGplot=BGplot;        %current's before nan
pct.DAproj=DAproj;          %current's before nan
pct.DAiso=ITplot;           %current's after  nan
%pct.rawIT=isub;
pct.Qplot=Q;
pct.Iox=nanmean(isub(oxid-2:oxid+2,:),1);

end
