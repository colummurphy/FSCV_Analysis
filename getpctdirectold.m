function pct = getpctdirect(isub, varargin)
%get pca directly from isub, and already found glitches (input argument)
smoothwidth=4;        %# of samples +/- to interpolate with nan around glitches/artifacts
argnum=1;
qthres=[];
ats=[];
cts=[];
includem=0;
rthresm=0.85;
rthresp=0.95;    %rthresout for ph
rthresb=0.95;  %rthresout for bg
imaxbg=[];
glitchids=[];
numpcs=3;       %da,ph,bg pc components to use default

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
        case 'includem'
            %include M (K(:,4)) in calculating dopamine or not
            includem=1;
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
  
KscaledM=[];        %movement templates repeated to find correlated data (ie. artifacts)
Mplot=[];           %also calculate M current component if we were to include all 4 signal categories
%if do not want to use movement component, but exists in loaded pca
%templates, remove from templates cts and ats variables
%still calculated KscaledM for artifact removal
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
end
dsdata=[];
if numpcs==3
    [Vc, F, Aproj, Qa]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
    E=isub-(Vc*Dproj);     %residuals
    Q=diag(E'*E)';      %nA^2
    IDAproj=Cu(1,:);            %current da projection
    pHplot=Cu(2,:);
    BGplot=Cu(3,:);
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    KDA=K(:,1);
    KPH=K(:,2);
    KBG=K(:,3);
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



%get correlation coefficients between data all time and movementtemplate
if exist('smoothdata')
    %if is function
    dsdata=smoothdata(isub,2,'gaussian',smoothwidth); 
else
    %use <2013 version of smooth data, need to check more carefully
    dsdata=smoothts(isub,'g',smoothwidth); 
end

rr=abs(corr(dsdata,KscaledM(:,1)));
rm=find(rr>rthresm);
rb=abs(corr(dsdata,KBG));
routb=find(rb>rthresb);
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
pct.pHplot=pHplot;         %current's before nan
pct.Mplot=Mplot;          %current's before nan
pct.BGplot=BGplot;        %current's before nan
pct.DAproj=DAproj;          %current's before nan
pct.DAiso=ITplot;           %current's after  nan
%pct.rawIT=isub;
pct.Qplot=Q;

end
