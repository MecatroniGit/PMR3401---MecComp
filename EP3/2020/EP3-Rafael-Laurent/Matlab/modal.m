function [A,w1,CC]=modal(nos_elementos,elementos,Kglobs,Mglobs,Lele)

%Função que realiza a análise modal

tic
CC=3*find(nos_elementos(:,2)==0 & (nos_elementos(:,1)==0 | nos_elementos(:,1)==4));

MCC=Mglobs;
KCC=Kglobs;
KCC(CC(2)-1,:)=[];
KCC(:,CC(2)-1)=[];
KCC(CC(1)-2:CC(1)-1,:)=[];
KCC(:,CC(1)-2:CC(1)-1)=[];

MCC(CC(2)-1,:)=[];
MCC(:,CC(2)-1)=[];
MCC(CC(1)-2:CC(1)-1,:)=[];
MCC(:,CC(1)-2:CC(1)-1)=[];

[A,V]=eigs(inv(MCC)*KCC,6,'smallestabs');

w1=sqrt(abs(V))/(2*pi);

toc
plotavibsmodes(A,nos_elementos,elementos,Lele)

end

