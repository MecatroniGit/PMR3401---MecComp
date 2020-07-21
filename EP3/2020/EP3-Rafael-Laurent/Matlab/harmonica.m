function harmonica(nos_elementos,MCC,KCC,nos_F1,nos_F2)

%Função que realiza a análise harmônica

modo=2;
Forcas=geraF(nos_F1,nos_F2,0,length(MCC),modo);

w=1:0.1:35;

noB=find(nos_elementos(:,1)==1.75 & nos_elementos(:,2)==4.5);
noC=find(nos_elementos(:,1)==1.75 & nos_elementos(:,2)==1.5);
vetB=zeros(length(w),2);
vetC=zeros(length(w),2);

for i=1:length(w)
    freq=(2*pi*w(i))^2;
    Mbarra=KCC-freq*MCC;
    solucao=Mbarra\Forcas;
    vetB(i,:)=[solucao(3*noB-2) solucao(3*noB-1)];
    vetC(i,:)=[solucao(3*noC-2) solucao(3*noC-1)];
end

modB=sqrt(vetB(:,1).^2+vetB(:,2).^2);
modC=sqrt(vetC(:,1).^2+vetC(:,2).^2);

figure
plot(w,modB)
hold on
plot(w,modC)
legend('Node B','Node C')
xlabel('Frequency (Hz)')
ylabel('Absolute Displacement (m)')
title('Frequency influence in Nodal Absolute Displacement')
saveas(gcf,'harmonicod.jpg')

end

