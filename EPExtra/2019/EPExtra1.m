clear all;
clc

%Criação da malha. Nós e elementos segundo a geometria do problema
disp("Creating Mesh");
h = 2;    %Tamanho do elemento triangular
r = 26;   %Raio do domínio

x = 0;
y = 0;
elements = [];

%Criação dos elementos triangulares a partir de coordenadas
while  y < r
    x=0;
    maxXl = sqrt(r^2 - y^2);    
    while x<= maxXl
        %Triangulo de baixo, na região do ar
        element.p1x = x;
        element.p1y = y;
        element.p2x = x+h;
        element.p2y = y;
        element.p3x = x+h;
        element.p3y = y+h;        
        elements = [elements element];
        
        %Triangulo de cima, na região do ar
        element.p1x = x;
        element.p1y = y;
        element.p2x = x+h;
        element.p2y = y+h;
        element.p3x = x;
        element.p3y = y+h;
        elements = [elements element];
        
        %Triangulo de cima, na região do solo
        element.p1x = x;
        element.p1y = -y;
        element.p2x = x+h;
        element.p2y = -y-h;
        element.p3x = x+h;
        element.p3y = -y;        
        elements = [elements element];
        
        %Triangulo de baixo, na região do solo
        element.p1x = x;
        element.p1y = -y;
        element.p2x = x;
        element.p2y = -y-h;
        element.p3x = x+h;
        element.p3y = -y-h;
        elements = [elements element];
        x = x+h; 
    end
    y = y+h;
end

node.x = 0;
node.y = 0;
nodes=[node];
%Definição do nós: Relaciona as coordenadas dos pontos dos triangulos
%Verifica a existência do ponto na lista. Se não encontrar, um novo nó é
%criado
for i = 1:length(elements)
    elements(i).xm =  (elements(i).p1x + elements(i).p2x + elements(i).p3x)/3;
    elements(i).ym =  (elements(i).p1y + elements(i).p2y + elements(i).p3y)/3;
    element = elements(i);
    %Vértice 1 do elemento
    elements(i).p1Node = -1;
    eNode.x = element.p1x;
    eNode.y = element.p1y;
    for j = 1:length(nodes) 
        if isequal(eNode, nodes(j))
            elements(i).p1Node = j;
            break
        end        
    end
    if elements(i).p1Node == -1
            nodes = [nodes, eNode];
            elements(i).p1Node = j+1;
    end
    %Vértice 2 do elemento
    elements(i).p2Node = -1;
    eNode.x = element.p2x;
    eNode.y = element.p2y;
    for j = 1:length(nodes)        
        if isequal(eNode, nodes(j))
            elements(i).p2Node = j;
        end
        
    end
    if elements(i).p2Node == -1
            nodes = [nodes, eNode];
            elements(i).p2Node = j+1;
    end
    %Vértice 3 do elemento
    elements(i).p3Node = -1;
    eNode.x = element.p3x;
    eNode.y = element.p3y;
    for j = 1:length(nodes)        
        if isequal(eNode, nodes(j))
            elements(i).p3Node = j;
        end        
    end
    if elements(i).p3Node == -1
            nodes = [nodes, eNode];
            elements(i).p3Node = j+1;
    end
end

%Aplica condições de contorno aos nós
disp("Apllying boundary conditions");
%Aplica condições de contorno e outros parâmetros nos elementos
border_threshold = 1;

%parâmetros dependentes do meio
miAr = 1.2566E-6;
miSolo = 2 * 1.2566E-6;
sigmaAr = 1E-10;
sigmaSolo = 1E-2;

Imax = 200;
rc = 0.02;

%Atribui potencial 0 nos pontos do contorno
for i = 1:length(nodes)
    nRad = sqrt((nodes(i).x)^2 + (nodes(i).y)^2);
    if nRad > (r - border_threshold)
        nodes(i).V = 0;
        nodes(i).A = 0;
        %disp("Added Contour");
    end
end

%Nós onde  estão as fontes (fios)
Vmax = 500;
Asource = - miAr * Imax /2/pi/rc;

bcNode1.x = 4;
bcNode1.y = 14;
bcNode2.x = 6;
bcNode2.y = 10;

for i = 1:length(nodes)
    %Procura pontos para as fontes
    if (nodes(i).x == bcNode1.x && nodes(i).y == bcNode1.y)||(nodes(i).x == bcNode2.x && nodes(i).y == bcNode2.y)
        nodes(i).V = Vmax;
        nodes(i).ASource = Asource;
        %disp("Added V source");
    end
    
    %Procura região do carro, atribuindo o valor do potencial
    if (nodes(i).x <=2 && nodes(i).y <=2 && nodes(i).y >=0)
        nodes(i).V = 0;
        %disp("Added car Surface");
    end
end


%Atribuição dos valores de mi e sigma aos elementos
for i = 1:length(elements)
    if elements(i).p3y > 0
        %disp("Element is air");
        elements(i).sigma = sigmaAr;
        elements(i).mi = miAr;
    else        
        elements(i).sigma = sigmaSolo;
        elements(i).mi = miSolo;
        %disp(elements(i));
    end
end
        


%Tratamento para o potencial Elétrico
%Construção da matriz global K
disp("Solving for V");
F = zeros(length(nodes),1);
K = zeros(length(nodes));
for i = 1:length(elements)    
    localMatrix = Ke(elements(i));  
    %disp(localMatrix);
    m = elements(i).p1Node;
    n = elements(i).p2Node;
    o = elements(i).p3Node;
    K(m,m) = K(m,m) + localMatrix(1,1);
    K(m,n) = K(m,n) + localMatrix(1,2);
    K(m,o) = K(m,o) + localMatrix(1,3);
    K(n,m) = K(n,m) + localMatrix(2,1);
    K(n,n) = K(n,n) + localMatrix(2,2);
    K(n,o) = K(n,o) + localMatrix(2,3);
    K(o,m) = K(o,m) + localMatrix(3,1);
    K(o,n) = K(o,n) + localMatrix(3,2);
    K(o,o) = K(o,o) + localMatrix(3,3);
end

%Contrução da matriz F, de carregamentos
for i = 1:length(nodes)
    if  isempty(nodes(i).V)
    else
        for j = 1:length(K)
            F(j) = F(j) - K(j,i) * nodes(i).V;
            K(i,j) = 0;
            K(j,i) = 0;                       
        end
        F(i) = nodes(i).V;
        K(i,i) = 1;       
    end
end

%Solução do sistema
V = linsolve(K, F);
%Atribuição dos valores calculados aos nós
for i = 1:length(nodes)
    nodes(i).V = V(i);
end

%Calculo das derivadas parciais
for i = 1:length(elements)
    [dx_,dy_] = Edxdy(elements(i), nodes);
    elements(i).Ex = -dx_;
    elements(i).Ey = -dy_;
end
    


%Exibição do resultado

VPlot = zeros(56/h, 30/h);
for k = 1:length(nodes)
  x = nodes(k).x;
  y = nodes(k).y;
  j = 1 + x/h;
  i = 1+(y+28)/h;
  VPlot(i,j) = V(k);
end
if h ==2 
    [X,Y] = meshgrid(0:h:28, -28:h:27); %h = 2
elseif h ==1
    [X,Y] = meshgrid(0:h:29, -28:h:28-h); %h = 1
end
X_ = [-fliplr(X) X];
Y_ = [fliplr(Y) Y];
figure;
VPlot_ = [fliplr(VPlot) VPlot];
mesh(X_, Y_, VPlot_);            
title("Potencial elétrico V");
xlabel("x [m]");
ylabel("y [m]");
zlabel("V [kV]");
figure;
contour(X_, Y_, VPlot_, 'ShowText','off', 'LabelSpacing',500);            
title("Potencial elétrico em curvas de nível");
c = colorbar;
c.Label.String = 'E [kV]'; %%Falta unidade
xlabel("x [m]");
ylabel("y [m]");
axis equal;
ylim([-5 28]);
hold on;
xms = zeros(length(elements),1);
yms = zeros(length(elements),1);
Exs = zeros(length(elements),1);
Eys = zeros(length(elements),1);
for i = 1:length(elements)
    xms(i) = elements(i).xm;
    yms(i) = elements(i).ym; 
    Exs(i) = elements(i).Ex;
    Eys(i) = elements(i).Ey;      
end
quiver([-fliplr(xms) xms],[fliplr(yms) yms],[-fliplr(Exs) Exs],[fliplr(Eys) Eys], 0.1);
legend({'V', "Exy"},'Location','southwest');

%Tratamento para o potencial vetor magnético
%Construção da matriz global KA
disp("Solving for A");
FA = zeros(length(nodes),1);
KA = zeros(length(nodes));
for i = 1:length(elements)    
    localMatrix = Ka(elements(i));  
    %disp(localMatrix);
    m = elements(i).p1Node;
    n = elements(i).p2Node;
    o = elements(i).p3Node;
    KA(m,m) = KA(m,m) + localMatrix(1,1);
    KA(m,n) = KA(m,n) + localMatrix(1,2);
    KA(m,o) = KA(m,o) + localMatrix(1,3);
    KA(n,m) = KA(n,m) + localMatrix(2,1);
    KA(n,n) = KA(n,n) + localMatrix(2,2);
    KA(n,o) = KA(n,o) + localMatrix(2,3);
    KA(o,m) = KA(o,m) + localMatrix(3,1);
    KA(o,n) = KA(o,n) + localMatrix(3,2);
    KA(o,o) = KA(o,o) + localMatrix(3,3);
end

%Contrução da matriz F, de carregamentos
for i = 1:length(nodes)
    %Fontes J
    if  isempty(nodes(i).ASource)
    else        
        FA(i) = nodes(i).ASource;             
    end
    %Borda (Dirichilet)
    if  isempty(nodes(i).A)
    else
        for j = 1:length(K)
            FA(j) = FA(j) - KA(j,i) * nodes(i).A;
            KA(i,j) = 0;
            KA(j,i) = 0;                       
        end
        FA(i) = nodes(i).V;
        KA(i,i) = 1;       
    end
end

%Solução do sistema
A = linsolve(KA, FA);
%Atribuição dos valores calculados aos nós
for i = 1:length(nodes)
    nodes(i).A = A(i);
end

%Calculo das derivadas:
%Densidade de fluxo magnético B
for i = 1:length(elements)
    [dx_,dy_] = Adxdy(elements(i), nodes);
    elements(i).Bx = dy_;
    elements(i).By = -dx_;
    elements(i).Hx = elements(i).Bx/elements(i).mi;
    elements(i).Hy = elements(i).By/elements(i).mi;
end

%Exibição do resultado
APlot = zeros(56/h, 30/h);
for k = 1:length(nodes)
  x = nodes(k).x;
  y = nodes(k).y;
  j = 1 + x/h;
  i = 1+(y+28)/h;
  APlot(i,j) = A(k);
end
if h ==2 
    [X,Y] = meshgrid(0:h:28, -28:h:27); %h = 2
elseif h ==1
    [X,Y] = meshgrid(0:h:29, -28:h:28-h); %h = 1
end
X_ = [-fliplr(X) X];
Y_ = [fliplr(Y) Y];
figure;
APlot_ = [fliplr(APlot) APlot];
mesh(X_, Y_, APlot_);            
title("Potencial Vetor Magnético A");
xlabel("x [m]");
ylabel("y [m]");
zlabel("A [Wb/m]");          %%Falta unidade
figure;
contour(X_, Y_, APlot_,16, 'ShowText','off', 'LabelSpacing',500);            
title("Potencial Vetor Magnético e Fluxo magnético");
c = colorbar;
c.Label.String = 'A [Wb/m]'; %%Falta unidade
xlabel("x [m]");
ylabel("y [m]");
axis equal;

Bxs = zeros(length(elements),1);
Bys = zeros(length(elements),1);
Hxs = zeros(length(elements),1);
Hys = zeros(length(elements),1);
for i = 1:length(elements)
    Bxs(i) = elements(i).Bx;
    Bys(i) = elements(i).By;    
    Hxs(i) = elements(i).Hx;
    Hys(i) = elements(i).Hy;  
end
hold on;
quiver([-fliplr(xms) xms],[fliplr(yms) yms],[fliplr(Bxs) Bxs],[-fliplr(Bys) Bys], 0.1);
legend({'A', "B [T]"},'Location','southwest');

figure;
quiver([-fliplr(xms) xms],[fliplr(yms) yms],[fliplr(Hxs) Hxs],[-fliplr(Hys) Hys], 0.1);
title("Campo magnético H ");
xlabel("x [m]");
ylabel("y [m]");
legend({'[Wb]'},'Location','southwest');
axis equal;


function [lMatrix] = Ke(element)
%Calculo do potencial eletrico para um elemento tringular
% 
x1 = element.p1x;
x2 = element.p2x;
x3 = element.p3x;
y1 = element.p1y;
y2 = element.p2y;
y3 = element.p3y;

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



lMatrix = [b1*b1+c1*c1 b1*b2+c1*c2 b1*b3+c1*c3;
           b1*b2+c1*c2 b2*b2+c2*c2 b2*b3+c2*c3;
           b1*b3+c1*c3 b2*b3+c2*c3 b3*b3+c3*c3];
       
lMatrix = 1/(4*Ae) * lMatrix;
lMatrix = (element.sigma)* lMatrix;
end

function [lMatrix] = Ka(element)
%Calculo do potencial eletrico para um elemento tringular
% 
x1 = element.p1x;
x2 = element.p2x;
x3 = element.p3x;
y1 = element.p1y;
y2 = element.p2y;
y3 = element.p3y;

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



lMatrix = [b1*b1+c1*c1 b1*b2+c1*c2 b1*b3+c1*c3;
           b1*b2+c1*c2 b2*b2+c2*c2 b2*b3+c2*c3;
           b1*b3+c1*c3 b2*b3+c2*c3 b3*b3+c3*c3];
       
lMatrix = 1/(4*Ae) * lMatrix;
lMatrix = element.mi * lMatrix;
end

function [dx,dy] = Edxdy(element, nodes)
%Calcula as derivadas parciais em um elemento
%   Detailed explanation goes here
x1 = element.p1x;
x2 = element.p2x;
x3 = element.p3x;
y1 = element.p1y;
y2 = element.p2y;
y3 = element.p3y;
V1 = nodes(element.p1Node).V;
V2 = nodes(element.p2Node).V;
V3 = nodes(element.p3Node).V;

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
end

function [dx,dy] = Adxdy(element, nodes)
%Calcula as derivadas parciais em um elemento
%   Detailed explanation goes here
x1 = element.p1x;
x2 = element.p2x;
x3 = element.p3x;
y1 = element.p1y;
y2 = element.p2y;
y3 = element.p3y;
A1 = nodes(element.p1Node).A;
A2 = nodes(element.p2Node).A;
A3 = nodes(element.p3Node).A;

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
end

