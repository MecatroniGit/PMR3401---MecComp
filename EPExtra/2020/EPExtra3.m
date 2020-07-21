close all;
clc

%Criação da malha. Nós e elementos segundo a geometria do problema
h = 2;    %Tamanho do elemento triangular
r = 26;   %Raio do domínio

x = 0;
y = 0;
elements = [];
elements=[];

%parâmetros dependentes do meio
miAr = 1.2566E-6;
miSolo = 2 * 1.2566E-6;
sigmaAr = 1E-10;
sigmaSolo = 1E-2;


%Criação dos elementos triangulares a partir de coordenadas
while  y < r
    x=0;
    maxXl = sqrt(r^2 - y^2);    
    while x<= maxXl
        
        %Triangulo de baixo, na região do ar        
        vetorajuda=zeros(1,13);
        vetorajuda(1:6)=[x y x+h y x+h y+h];
        elements=[elements;vetorajuda];

        %Triangulo de cima, na região do ar        
        vetorajuda=zeros(1,13);
        vetorajuda(1:6)=[x y x+h y+h x y+h];
        elements=[elements;vetorajuda];
        

        %Triangulo de cima, na região do solo
        vetorajuda=zeros(1,13);
        vetorajuda(1:6)=[x -y x+h -y-h x+h -y];
        elements=[elements;vetorajuda];
     

        %Triangulo de baixo, na região do solo        
        vetorajuda=zeros(1,13);
        vetorajuda(1:6)=[x -y x -y-h x+h -y-h];
        elements=[elements;vetorajuda];
        
        x = x+h; 
    end
    y = y+h;
end

nos=[];
%Definição do nós: Relaciona as coordenadas dos pontos dos triangulos
%Verifica a existência do ponto na lista. Se não encontrar, um novo nó é
%criado
for i = 1:length(elements)
    %Calcula o xm
    elements(i,7)=(elements(i,1)+elements(i,3)+...
        elements(i,5))/3;
    %Calcula o ym
    elements(i,8)=(elements(i,2)+elements(i,4)+...
        elements(i,6))/3;
    
    
    %Vértice 1 do elemento    
    if length(nos)>0
        p1ehno=find(nos(:,1)==elements(i,1) &...
            nos(:,2)==elements(i,2));
    else
        p1ehno=[];
    end
    if length(p1ehno)==0
        nos(end+1,:)=[elements(i,1)  elements(i,2)];
        %Miguezinho que eu tive que dar
        %Se conseguirem, mudem o migue
        %Da pra mudar esse if pra if length(nos)==2
        %Ou então if elements(1,9)==1
        %Ou até mesmo calcula p1ehno de novo e faz if p1ehno==1
        %Tudo funciona, mas é um migue, mudem os migues
        if i==1
            elements(i,9)=1;
        else
            elements(i,9)=length(nos);
        end
    else
        elements(i,9)=p1ehno;
    end

    %Vértice 2 do elemento
    p2ehno=find(nos(:,1)==elements(i,3) &...
        nos(:,2)==elements(i,4));
    if length(p2ehno)==0
        nos(end+1,:)=[elements(i,3)  elements(i,4)];
        elements(i,10)=length(nos);
    else
        elements(i,10)=p2ehno;
    end    
    
    %Vértice 3 do elemento    
    p3ehno=find(nos(:,1)==elements(i,5) &...
        nos(:,2)==elements(i,6));
    if length(p3ehno)==0
        nos(end+1,:)=[elements(i,5)  elements(i,6)];
        elements(i,11)=length(nos);
    else
        elements(i,11)=p3ehno;
    end
    %Atribuição dos valores de mi e sigma aos elementos
    if elements(i,6)>0
        elements(i,12)=sigmaAr;
        elements(i,13)=miAr;
    else
        elements(i,12)=sigmaSolo;
        elements(i,13)=miSolo;
    end
end


%Aplica condições de contorno e outros parâmetros nos elementos
limite = 1;

Imax = 200;
rc = 0.02;


%Atribui potencial 0 nos pontos do contorno
nos2=zeros(length(nos),2);
nos2(:,:)=-1;
nos=[nos nos2];
zero=find(nos(:,1).^2+nos(:,2).^2>(r-limite)^2);
nos(zero,3:4)=0;

%Nós onde  estão as fontes (fios)
Vmax = 500;
Asource = - miAr * Imax /2/pi/rc;

%Procura pontos para as fontes
Asource=find((nos(:,1)==6 & nos(:,2)==10) | nos(:,1)==4 & nos(:,2)==14);
nos(Asource,3)=Vmax;

%Procura região do carro, atribuindo o valor do potencial
Carro=find(nos(:,1)<=2 & nos(:,2)<=1.5 & nos(:,2)>=0);
nos(Carro,3)=0;

%Tratamento para o potencial Elétrico
%Construção da matriz global K
K=zeros(length(nos));
for i=1:length(elements)
    
    x1 = elements(i,1);
    y1 = elements(i,2);
    x2 = elements(i,3);
    y2 = elements(i,4);
    x3 = elements(i,5);
    y3 = elements(i,6);

    a1= x2*y3-x3*y2;
    a2 = x3*y1 - x1*y3;
    a3 = x1*y2 - x2*y1;

    b1= y2 - y3;
    b2= y3 - y1;
    b3= y1 - y2;

    c1= x3 - x2;
    c2= x1 - x3;
    c3= x2 - x1;
    Ae= (b1*c2 - b2*c1)/2;



    Kele = [b1*b1+c1*c1 b1*b2+c1*c2 b1*b3+c1*c3;
               b1*b2+c1*c2 b2*b2+c2*c2 b2*b3+c2*c3;
               b1*b3+c1*c3 b2*b3+c2*c3 b3*b3+c3*c3];

    Kele = 1/(4*Ae) * Kele;
    Kele = elements(i,12)* Kele;
    m = elements(i,9);
    n = elements(i,10);
    o = elements(i,11);
    K(m,m) = K(m,m) + Kele(1,1);
    K(m,n) = K(m,n) + Kele(1,2);
    K(m,o) = K(m,o) + Kele(1,3);
    K(n,m) = K(n,m) + Kele(2,1);
    K(n,n) = K(n,n) + Kele(2,2);
    K(n,o) = K(n,o) + Kele(2,3);
    K(o,m) = K(o,m) + Kele(3,1);
    K(o,n) = K(o,n) + Kele(3,2);
    K(o,o) = K(o,o) + Kele(3,3);
end

F=zeros(length(nos),1);
%Contrução da matriz F, de carregamentos
for i = 1:length(nos)
    if  nos(i,3)~=-1
        for j = 1:length(K)
            F(j) = F(j) - K(j,i) * nos(i,3);
            K(i,j) = 0;
            K(j,i) = 0;                       
        end
        F(i) = nos(i,3);
        K(i,i) = 1;       
    end
end

%Solução do sistema
V = linsolve(K, F);
%Atribuição dos valores calculados aos nós
for i = 1:length(nos)
    nos(i,3) = V(i);
end

%Calculo das derivadas parciais
E=zeros(length(elements),2);
for i = 1:length(elements)
    x1 = elements(i,1);
    y1 = elements(i,2);
    
    x2 = elements(i,3);
    y2 = elements(i,4);
    
    x3 = elements(i,5);
    y3 = elements(i,6);
    
    V1 = nos(elements(i,9),3);
    V2 = nos(elements(i,10),3);
    V3 = nos(elements(i,11),3);

    a1= x2*y3-x3*y2;
    a2 = x3*y1 - x1*y3;
    a3 = x1*y2 - x2*y1;

    b1= y2 - y3;
    b2= y3 - y1;
    b3= y1 - y2;

    c1= x3 - x2;
    c2= x1 - x3;
    c3= x2 - x1;
    Ae= (b1*c2 - b2*c1)/2;


    dx = (1/(2*Ae))* (b1*V1 + b2*V2 + b3*V3);
    dy = (1/(2*Ae))* (c1*V1 + c2*V2 + c3*V3);
    E(i,1)= -dx;
    E(i,2)= -dy;
end
    
%Exibição do resultado
VPlot = zeros(56/h, 30/h);
for k = 1:length(nos)
  x = nos(k,1);
  y = nos(k,2);
  j = 1 + x/h;
  i = 1+(y+28)/h;
  VPlot(i,j) = V(k);
end
if h ==2 
    [X,Y] = meshgrid(-28:h:28, -28:h:27); %h = 2
elseif h ==1
    [X,Y] = meshgrid(-29:h:29, -28:h:28-h); %h = 1
end

figure;
VPlot_ = [fliplr(VPlot) VPlot];
VPlot_(:,length(VPlot_)/2)=[];

surf(X, Y, VPlot_);            
title("Potencial elétrico  2 V");
xlabel("x [m]");
ylabel("y [m]");
zlabel("V [kV]");

figure;
contour(X, Y, VPlot_, 'ShowText','off', 'LabelSpacing',500);            
title("Potencial elétrico em curvas de nível");
c = colorbar;
c.Label.String = 'E [kV]';
xlabel("x [m]");
ylabel("y [m]");

axis equal;
ylim([-5 28]);
hold on;

quiver([-fliplr(elements(:,7)') elements(:,7)'],...
[fliplr(elements(:,8)') elements(:,8)'],...
[-fliplr(E(:,1)') E(:,1)'],[fliplr(E(:,2)') E(:,2)'],2);
legend({'V ', "Exy "},'Location','southwest');

%Tratamento para o potencial vetor magnético
%Construção da matriz global KA
KA = zeros(length(nos));
for i = 1:length(elements)    
    x1 = elements(i,1);
    x2 = elements(i,3);
    x3 = elements(i,5);
    y1 = elements(i,2);
    y2 = elements(i,4);
    y3 = elements(i,6);

    a1= x2*y3-x3*y2;
    a2 = x3*y1 - x1*y3;
    a3 = x1*y2 - x2*y1;

    b1= y2 - y3;
    b2= y3 - y1;
    b3= y1 - y2;

    c1= x3 - x2;
    c2= x1 - x3;
    c3= x2 - x1;
    Ae= (b1*c2 - b2*c1)/2;

    Kele = [b1*b1+c1*c1 b1*b2+c1*c2 b1*b3+c1*c3;
               b1*b2+c1*c2 b2*b2+c2*c2 b2*b3+c2*c3;
               b1*b3+c1*c3 b2*b3+c2*c3 b3*b3+c3*c3];

    Kele = 1/(4*Ae) * Kele;
    Kele = elements(i,13) * Kele;    %disp(localMatrix);
    m = elements(i,9);
    n = elements(i,10);
    o = elements(i,11);
    KA(m,m) = KA(m,m) + Kele(1,1);
    KA(m,n) = KA(m,n) + Kele(1,2);
    KA(m,o) = KA(m,o) + Kele(1,3);
    KA(n,m) = KA(n,m) + Kele(2,1);
    KA(n,n) = KA(n,n) + Kele(2,2);
    KA(n,o) = KA(n,o) + Kele(2,3);
    KA(o,m) = KA(o,m) + Kele(3,1);
    KA(o,n) = KA(o,n) + Kele(3,2);
    KA(o,o) = KA(o,o) + Kele(3,3);
end

%Contrução da matriz F, de carregamentos
FA = zeros(length(nos),1);
for i = 1:length(nos)
    %Fontes J
    if  length(find(Asource==i))~=0
        FA(i) = - miAr * Imax /2/pi/rc;             
    end
    %Borda (Dirichilet)
    if  nos(i,4)~=-1
        for j = 1:length(K)
            FA(j) = FA(j) - KA(j,i) * nos(i,4);
            KA(i,j) = 0;
            KA(j,i) = 0;                       
        end
        FA(i) = nos(i,3);
        KA(i,i) = 1;       
    end
end


%Solução do sistema
A = linsolve(KA, FA);
%Atribuição dos valores calculados aos nós
for i = 1:length(nos)
    nos(i,4) = A(i);
end

%Calculo das derivadas:
%Densidade de fluxo magnético B e campo magnético H
B=zeros(length(elements),2);
H=zeros(length(elements),2);
for i = 1:length(elements)
    x1 = elements(i,1);
    x2 = elements(i,3);
    x3 = elements(i,5);
    y1 = elements(i,2);
    y2 = elements(i,4);
    y3 = elements(i,6);
    
    A1 = nos(elements(i,9),4);
    A2 = nos(elements(i,10),4);
    A3 = nos(elements(i,11),4);

    a1= x2*y3-x3*y2;
    a2 = x3*y1 - x1*y3;
    a3 = x1*y2 - x2*y1;

    b1= y2 - y3;
    b2= y3 - y1;
    b3= y1 - y2;

    c1= x3 - x2;
    c2= x1 - x3;
    c3= x2 - x1;
    Ae= (b1*c2 - b2*c1)/2;

    dx = (1/(2*Ae))* (b1*A1 + b2*A2 + b3*A3);
    dy = (1/(2*Ae))* (c1*A1 + c2*A2 + c3*A3);
    B(i,1)=dy;
    B(i,2)=-dx;
    
    H(i,1)=B(i,1)/elements(i,13);
    H(i,2)=B(i,2)/elements(i,13);
end

%Exibição do resultado
APlot = zeros(56/h, 30/h);
for k = 1:length(nos)
  x = nos(k,1);
  y = nos(k,2);
  j = 1 + x/h;
  i = 1+(y+28)/h;
  APlot(i,j) = A(k);
end
if h ==2 
    [X,Y] = meshgrid(-28:h:28, -28:h:27); %h = 2
elseif h ==1
    [X,Y] = meshgrid(-29:h:29, -28:h:28-h); %h = 1
end

APlot_ = [fliplr(APlot) APlot];
APlot_(:,length(APlot_)/2)=[];     

figure;
surf(X, Y, APlot_);            
title("Potencial Vetor Magnético A ");
xlabel("x [m]");
ylabel("y [m]");
zlabel("A [Wb/m]");   


figure;
contour(X, Y, APlot_,16, 'ShowText','off', 'LabelSpacing',500);            
title("Potencial Vetor Magnético e Fluxo magnético ");
c = colorbar;
c.Label.String = 'A [Wb/m]'; %%Falta unidade
xlabel("x [m]");
ylabel("y [m]");
axis equal;

hold on;
quiver([-fliplr(elements(:,7)') elements(:,7)'],...
[fliplr(elements(:,8)') elements(:,8)'],...
[fliplr(B(:,1)') B(:,1)'],[-fliplr(B(:,2)') B(:,2)'], 2);
legend({'A ', "B [T] "},'Location','southwest');

figure;
quiver([-fliplr(elements(:,7)') elements(:,7)'],...
[fliplr(elements(:,8)') elements(:,8)'],...
[fliplr(H(:,1)') H(:,1)'],[-fliplr(H(:,2)') H(:,2)'], 2);
title("Campo magnético H ");
xlabel("x [m]");
ylabel("y [m]");
legend({'[Wb]'},'Location','southwest');
axis equal;


