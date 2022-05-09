function [Vc F Aproj Qa] = generatePCR(Ats,Cts)
dir=[cd '\data\'];
%Ats = training set voltammetric matrix nxm, n = # potential steps, m = #
%training set samples
%Cts = jxm training set reference conc values (j = # analytes, m = #
%training set samples)
dir=cd;
%cd('C:\Users\Dai5H\Documents\MATLAB\FCV\')
%load('Concentration matrix MIT.txt') % 2x14 ??? Weights?
%load('CV matrix MIT.txt') % 175x14, first 175x3 corresponds to 3 conc of DA? Training set? 8-13 columns are pH, first 7 DA
CV_matrix=Ats;

meanmat=repmat(mean(CV_matrix'),size(CV_matrix,2),1)'; % 1x175 array representing mean of entire CV_matrix tiled 14 times

A=CV_matrix-meanmat;            %Start using this 08/11/2017, stopped using when using Patra templates 03/30/2018
%A=CV_matrix;    %nxm
[U S V]=svd(A); %singular value decomposition
k=4;        %# of relevant PC's, for BG/Movement
k=size(Cts,1);
Vc=U(:,1:k);        %data matrix contiaining all relevant PC's
%Vc=U(:,1:end);        %data matrix contiaining all relevant PC's

Aproj=Vc'*A;    %projections of training set data for current at each voltage step onto relevant PC's
F=Cts*Aproj'*inv(Aproj*Aproj'); %regression matrix relating data projections to concentrations

Qa=calculateQa(S);

end