clear all;
close all;
clc

% Constantes
h = 2;    
r = 26; 
miAr = 1.2566E-6;
miSolo = 2 * 1.2567E-6;
sigmaAr = 1E-10;
sigmaSolo = 1E-2;
Imax = 200;
rc = 0.02;
Vmax = 500;

% Nós onde  estão as fontes
Asource = - miAr * Imax/(2* pi*rc);
bcNode1.x = 4;
bcNode1.y = 14;
bcNode2.x = 6;
bcNode2.y = 10;

% Criação da malha
elements = elementCreation(r,h);
[nodes, elements] = nodeCreation(elements);

% Potencial 0 nos pontos da borda
for i = 1:length(nodes)
    node = nodes(i);
    nRad = sqrt((node.x)^2 + (node.y)^2);
    if nRad >= r
        nodes(i).V = 0;
        nodes(i).A = 0;
    end
    % Procura pontos para as fontes
    if (nodes(i).x == bcNode1.x && nodes(i).y == bcNode1.y)||(nodes(i).x == bcNode2.x && nodes(i).y == bcNode2.y)
        nodes(i).V = Vmax;
        nodes(i).ASource = Asource;
    end 
    % Procura região do carro
    if (nodes(i).x <=2 && nodes(i).y <=1.5 && nodes(i).y >=0)
        nodes(i).V = 0;
    end
end

%Atribuição dos valores de mi e sigma aos elementos
for i = 1:length(elements)
    if elements(i).p3y > 0
        elements(i).sigma = sigmaAr;
        elements(i).mi = miAr;
    else        
        elements(i).sigma = sigmaSolo;
        elements(i).mi = miSolo;
    end
end

% Potencial Elétrico e Magnético
F = zeros(length(nodes),1);
K = globalMatrix(elements,length(nodes),1);
FA = zeros(length(nodes),1);
KA = globalMatrix(elements,length(nodes),0);

% Matrizes F
for i = 1:length(nodes)
    % Matriz FV
    if  ~isempty(nodes(i).V)
        for j = 1:length(K)
            F(j) = F(j) - K(j,i) * nodes(i).V;
            K(i,j) = 0;
            K(j,i) = 0;                       
        end
        F(i) = nodes(i).V;
        K(i,i) = 1;       
    end
    % Matriz FA
    if  ~isempty(nodes(i).ASource)       
        FA(i) = nodes(i).ASource;             
    end
    % Borda
    if  ~isempty(nodes(i).A)
        for j = 1:length(K)
            FA(j) = FA(j) - KA(j,i) * nodes(i).A;
            KA(i,j) = 0;
            KA(j,i) = 0;                       
        end
        FA(i) = nodes(i).V;
        KA(i,i) = 1;       
    end
end

% Solução
V = linsolve(K, F);
A = linsolve(KA, FA);

% Atribuição dos valores calculados aos nós
for i = 1:length(nodes)
    nodes(i).V = V(i);
    nodes(i).A = A(i);
end

% Cálculo das derivadas
for i = 1:length(elements)
    [dx_,dy_] = Edxdy(elements(i), nodes, 1);
    elements(i).Ex = -dx_;
    elements(i).Ey = -dy_;
    [dxA_,dyA_] = Edxdy(elements(i), nodes, 0);
    elements(i).Bx = dyA_;
    elements(i).By = -dxA_;
    elements(i).Hx = elements(i).Bx/elements(i).mi;
    elements(i).Hy = elements(i).By/elements(i).mi;
end

% Plots
VPlot = zeros(56/h, 30/h);
APlot = zeros(56/h, 30/h);
for k = 1:length(nodes)
  x = nodes(k).x;
  y = nodes(k).y;
  j = 1 + x/h;
  i = 1+(y+28)/h;
  VPlot(i,j) = V(k);
  APlot(i,j) = A(k);
end
VPlot_ = [fliplr(VPlot) VPlot];
APlot_ = [fliplr(APlot) APlot];

[X,Y] = meshgrid(0:h:28, -28:h:27); 
X_ = [-fliplr(X) X];
Y_ = [fliplr(Y) Y];

% V
figure;
contour(X_, Y_, VPlot_, 'ShowText','off', 'LabelSpacing',500);            
title("V");
c = colorbar;
c.Label.String = 'V [kV]';
xlabel("x [m]");
ylabel("y [m]");
axis equal;
ylim([-5 28]);
legend("V",'Location','southwest');

% E
figure();
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
title("E");
xlabel("x [m]");
ylabel("y [m]");
legend("E",'Location','southwest');

% A
figure;
contour(X_, Y_, APlot_,16, 'ShowText','off', 'LabelSpacing',500);            
title("A");
c = colorbar;
c.Label.String = 'A'; 
xlabel("x [m]");
ylabel("y [m]");
axis equal;
legend("A",'Location','southwest');

% B
figure();
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
quiver([-fliplr(xms) xms],[fliplr(yms) yms],[fliplr(Bxs) Bxs],[-fliplr(Bys) Bys], 0.1);
title("B");
xlabel("x [m]");
ylabel("y [m]");
legend("B",'Location','southwest');

% H
figure;
quiver([-fliplr(xms) xms],[fliplr(yms) yms],[fliplr(Hxs) Hxs],[-fliplr(Hys) Hys], 0.1);
title("H ");
xlabel("x [m]");
ylabel("y [m]");
legend('H','Location','southwest');
axis equal;

function [K_el] = Ke(element, parameter)
    % Montagem da matriz K do elemento
    x1 = element.p1x;
    x2 = element.p2x;
    x3 = element.p3x;
    y1 = element.p1y;
    y2 = element.p2y;
    y3 = element.p3y;

    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;
    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;
    
    Ae = (b1*c2 - b2*c1)/2;
    
    if parameter
        const = element.sigma;
    else
        const = element.mi;
    end 
    
    K_el = 1/(4*Ae)*const*[b1*b1+c1*c1 b1*b2+c1*c2 b1*b3+c1*c3;
                           b1*b2+c1*c2 b2*b2+c2*c2 b2*b3+c2*c3;
                           b1*b3+c1*c3 b2*b3+c2*c3 b3*b3+c3*c3];
end

function elements = elementCreation(r,h)
    x = 0;
    y = 0;
    elements = [];

    %Criação dos elementos triangulares a partir de coordenadas
    while  y < r
        x=0;
        maxXl = sqrt(r^2 - y^2);    
        while x<= maxXl
            %Triangulo de baixo, na região do ar
            element1.p1x = x;
            element1.p1y = y;
            element1.p2x = x+h;
            element1.p2y = y;
            element1.p3x = x+h;
            element1.p3y = y+h;       
            %Triangulo de cima, na região do ar
            element2.p1x = x;
            element2.p1y = y;
            element2.p2x = x+h;
            element2.p2y = y+h;
            element2.p3x = x;
            element2.p3y = y+h;
            %Triangulo de cima, na região do solo
            element3.p1x = x;
            element3.p1y = -y;
            element3.p2x = x+h;
            element3.p2y = -y-h;
            element3.p3x = x+h;
            element3.p3y = -y;        
            %Triangulo de baixo, na região do solo
            element4.p1x = x;
            element4.p1y = -y;
            element4.p2x = x;
            element4.p2y = -y-h;
            element4.p3x = x+h;
            element4.p3y = -y-h;
            
            elements = [elements element1 element2 element3 element4];
            x = x+h; 
        end
        y = y+h;
    end
end

function [nodes, elements] = nodeCreation(elements)
    node.x = 0;
    node.y = 0;
    nodes=[node];
    %Definição do nós: Relaciona as coordenadas dos pontos dos triangulos
    %Verifica a existência do ponto na lista. Se não encontrar, um novo nó é criado
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
    
end

function [dx,dy] = Edxdy(element, nodes, parameter)
%Calcula as derivadas parciais em um elemento
    x1 = element.p1x;
    x2 = element.p2x;
    x3 = element.p3x;
    y1 = element.p1y;
    y2 = element.p2y;
    y3 = element.p3y;
    
    if parameter
        V1 = nodes(element.p1Node).V;
        V2 = nodes(element.p2Node).V;
        V3 = nodes(element.p3Node).V;
    else
        V1 = nodes(element.p1Node).A;
        V2 = nodes(element.p2Node).A;
        V3 = nodes(element.p3Node).A;
    end
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

function K = globalMatrix(elements, lenNodes, parameter)
    K = zeros(lenNodes);
    for i = 1:length(elements)    
        localMatrix = Ke(elements(i),parameter);
        n1 = elements(i).p1Node;
        n2 = elements(i).p2Node;
        n3 = elements(i).p3Node;
        K(n1,n1) = K(n1,n1) + localMatrix(1,1);
        K(n1,n2) = K(n1,n2) + localMatrix(1,2);
        K(n1,n3) = K(n1,n3) + localMatrix(1,3);
        K(n2,n1) = K(n2,n1) + localMatrix(2,1);
        K(n2,n2) = K(n2,n2) + localMatrix(2,2);
        K(n2,n3) = K(n2,n3) + localMatrix(2,3);
        K(n3,n1) = K(n3,n1) + localMatrix(3,1);
        K(n3,n2) = K(n3,n2) + localMatrix(3,2);
        K(n3,n3) = K(n3,n3) + localMatrix(3,3);
    end 
end