function [MCC,KCC,nos_F1,nos_F2] = transiente(nos_elementos,CC,Mglobs,Kglobs)

%Função que realiza a análise transiente
global Elementos alfa beta gama beta2 L1 L2 L3 E

nos_F1=find(nos_elementos(:,2)==14);
nos_F1=nos_F1(1:Elementos:end);
nos_F2=find(nos_elementos(:,1)==1);
nos_F2=nos_F2(1:Elementos:end);

ti=0;
dt=0.02;
tf=10;

MCC=Mglobs;
KCC=Kglobs;


KCC(CC(2)-1,:)=0;
KCC(:,CC(2)-1)=0;
KCC(CC(2)-1,CC(2)-1)=1;

KCC(CC(1)-2:CC(1)-1,:)=0;
KCC(:,CC(1)-2:CC(1)-1)=0;
KCC(CC(1)-1,CC(1)-1)=1;
KCC(CC(1)-2,CC(1)-2)=1;

MCC(CC(2)-1,:)=0;
MCC(:,CC(2)-1)=0;
MCC(CC(2)-1,CC(2)-1)=1;

MCC(CC(1)-2:CC(1)-1,:)=0;
MCC(:,CC(1)-2:CC(1)-1)=0;
MCC(CC(1)-1,CC(1)-1)=1;
MCC(CC(1)-2,CC(1)-2)=1;

Cglobs=alfa*MCC+beta*KCC;

Forcas=zeros(length(Mglobs),length(ti:dt:tf));
desloc=zeros(length(Mglobs),length(ti:dt:tf));
vel=zeros(length(Mglobs),length(ti:dt:tf));
acel=zeros(length(Mglobs),length(ti:dt:tf));

modo=1; %transiente
Forcas(:,1)=geraF(nos_F1,nos_F2,0,length(Mglobs),modo);

acel(:,1)=inv(MCC)*(Forcas(:,1)-Cglobs*vel(:,1)-KCC*desloc(:,1));
Mbarra=MCC+dt*gama*Cglobs+dt^2*beta2*KCC;

for t=(ti+dt):dt:tf %A cada iteração calcula a{n+1=t} e F{n+1=t}
    Forcas(:,round(t/dt)+1)=geraF(nos_F1,nos_F2,t,length(MCC),modo);
    Mbarra=MCC+dt*gama*Cglobs+dt^2*beta2*KCC;
    
    Fbarra=Forcas(:,round(t/dt)+1)-Cglobs*(vel(:,round(t/dt))+...
    dt*(1-gama)*acel(:,round(t/dt)))-KCC*(desloc(:,round(t/dt))+...
    dt*vel(:,round(t/dt))+dt^2/2*(1-2*beta2)*acel(:,round(t/dt)));
    
    acel(:,round(t/dt)+1)=Mbarra\Fbarra;
    vel(:,round(t/dt)+1)=vel(:,round(t/dt))+...
    dt*((1-gama)*acel(:,round(t/dt))+gama*acel(:,round(t/dt)+1));
    desloc(:,round(t/dt)+1)=desloc(:,round(t/dt))+dt*vel(:,round(t/dt))...
    +dt^2/2*((1-2*beta2)*acel(:,round(t/dt))+2*beta2*acel(:,round(t/dt)+1));
    
    desloc(CC(2)-1,round(t/dt)+1)=0;
    desloc(CC(1)-2:CC(1)-1,round(t/dt)+1)=0;
    vel(CC(2)-1,round(t/dt)+1)=0;
    vel(CC(1)-2:CC(1)-1,round(t/dt)+1)=0;
    acel(CC(2)-1,round(t/dt)+1)=0;
    acel(CC(1)-2:CC(1)-1,round(t/dt)+1)=0;
end

noA=find(nos_elementos(:,1)==1 & nos_elementos(:,2)==14);
noB=find(nos_elementos(:,1)==1.75 & nos_elementos(:,2)==4.5);
noD=find(nos_elementos(:,1)==0.25 & nos_elementos(:,2)==1.5);
noE=find(nos_elementos(:,1)==3 & nos_elementos(:,2)==7);

%Plota os deslocamentos
figure()
plot(ti:dt:tf,desloc(3*noD-2,:)) %deslocamento em x
hold on
plot(ti:dt:tf,desloc(3*noE-2,:)) %deslocamento em x
legend('Node D','Node E')
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Nodes displacement in the X axis')
saveas(gcf,'transemx.jpg')

figure()
plot(ti:dt:tf,desloc(3*noD-1,:)) %deslocamento em y
hold on
plot(ti:dt:tf,desloc(3*noE-1,:)) %deslocamento em y
legend('Node D','Node E')
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Nodes displacement in the Y axis')
saveas(gcf,'transemy.jpg')

figure()
plot(ti:dt:tf,sqrt(desloc(3*noD-1,:).^2+desloc(3*noD-2,:).^2))
hold on
plot(ti:dt:tf,sqrt(desloc(3*noE-1,:).^2+desloc(3*noE-2,:).^2))
legend('Node D','Node E')
xlabel('Time (s)')
ylabel('Absolute Displacement (m)')
title('Nodes absolute displacement value')
saveas(gcf,'transemabs.jpg')

%Plota as tensões
%Tensão versão 1
dxA=E*(desloc(3*noA-2,:)-desloc(3*find(nos_elementos(:,1)==1+L1/...
    Elementos & nos_elementos(:,2)==14)-2,:))/(L1/Elementos);
tensaoA=abs(dxA);

dxB=E*(desloc(3*noB-2,:)-desloc(3*find(nos_elementos(:,1)==1.75+2.5/...
    Elementos & nos_elementos(:,2)==4.5+L2/Elementos)-2,:))/(sqrt(L2^2+(L3-1.5)^2)/Elementos);
dyB=E*(desloc(3*noB-1,:)-desloc(3*find(nos_elementos(:,1)==1.75-2.5/Elementos &...
    nos_elementos(:,2)==4.5-L2/Elementos)-1,:))/(sqrt(L2^2+(L3-1.5)^2)/Elementos);

tensaoB=sqrt(dxB.^2+dyB.^2);

figure()
plot(ti:dt:tf,tensaoA)
hold on
plot(ti:dt:tf,tensaoB)
xlabel('Time (s)')
ylabel('Stress (Pa)')
title('Stress Curve along the time')
legend('Node A','Node B')

end

