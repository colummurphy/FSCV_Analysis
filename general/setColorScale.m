%getColorPlotTemplates

function [map, cmin]=setColorScale(cmax)

%cmax=3;
cmin=-2/3*cmax;
%tick1=cmin/2.68;
tick1length=126;
tick2length=126;
ticklength=50;
%samplerate=sampling_freq;          %10samples/s

map=[0.4 0.1 0.5; 0.3 0.2 0.5; 0.1 0.2 0.3];
map=[0 0 0; 0 0 1; 0 0 0; 1 1 0; 1 0 1; 0 0.4 0; 0 1 0; 0 0 0];
cscale=(1:-(1/(ticklength-15)):0);
%cscaley=(.95:-((.95-0.1)/ticklength):0.1);
cscaley=(1:-((1-0.03)/ticklength):0.03);
%cscaley2=(1:-((1-0.25)/ticklength):0.25);
cscaley2=(1:-((1-0.35)/ticklength):0.35);
cscaley3=(0:((.25)/ticklength):0.25);
cscaleblue=(1:-(1/(ticklength+12)):0);
cblue2=[zeros(length(cscaleblue),1) zeros(length(cscaleblue),1) cscaleblue'];
%cyellow=[cscaley2' cscaley' zeros(length(cscale),1)];
cyellow=[cscaley2' cscaley' cscaley3'];
%ccool=flipud(cool(round(ticklength/3)));
ccool =[

    1.0000         0    1.0000;
    0.9375    0.0625    1.0000;
    0.8750    0.1250    1.0000;
    0.8125    0.1875    1.0000;
    0.7500    0.2500    1.0000;
    0.6875    0.3125    1.0000;
    0.6250    0.3750    1.0000;
    0.5625    0.4375    1.0000;
    0.5000    0.5000    1.0000;
    0.4375    0.5625    1.0000;
    0.3750    0.6250    1.0000;
    0.3125    0.6875    1.0000;
    0.2500    0.7500    1.0000;
    0.1875    0.8125    1.0000;
    0.1250    0.8750    1.0000;
    0.0625    0.9375    1.0000;
         0    1.0000    1.0000;
         ];
cscalehalf=(1:-(1/ticklength*1.5):0);
cscalehalf2=(1:-(0.5/ticklength*1.5):0.5);
btodgreen=[zeros(length(cscalehalf),1)  cscalehalf2' cscalehalf'];
cscalehalf3=(0.5:(1/ticklength*.5):1);
gtolime=[zeros(length(cscalehalf3),1) cscalehalf3' zeros(length(cscalehalf3),1)];
black=zeros(5,3);
red=[1 0 0];
map=[black; cblue2;  red; 0.65 0 0; cyellow; ccool; btodgreen; gtolime; black; black];
        %range=(cmin:double(1/length(map)):cmax);
        range=[double(cmin), double(.5*cmin), double(.25*cmin), 0, double(1/4*cmax), double(2/4*cmax), double(3/4*cmax), cmax]; 
        %range=[range,isarmax2];

end