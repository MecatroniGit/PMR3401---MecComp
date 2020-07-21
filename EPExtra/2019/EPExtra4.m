clear all
close all
clc

%% Parâmetros
global sigma_ar sigma_solo mi_ar mi_solo
mi_ar = 1.2566*10^-6;
mi_solo = 2*1.2566*10^-6;
sigma_ar = 10^-10;
sigma_solo = 10^-2;
h = input('passo da malha (1 ou 2): '); 

%% Construção do modelo de elementos finitos
[nodes, Nx] = CreateNodes(h);
[K] = CreateElements(h, nodes, Nx);
[KV, KB, FV, FB] = BoundaryConditions(h, K, nodes, Nx, mi_ar);

disp("Matrizes obtidas");

%% Resolução dos sistemas
Vnodes = linsolve(KV,FV); % Tensão nos nós
Aznodes = linsolve(KB,FB);% Potencial Magnético nos nós

disp("Sistema Resolvido");

%% Transforma os valores nodais em uma malha X, Y
V = zeros(52/h, 26/h);  % Matriz das tensões, em X e Y
Az = zeros(52/h, 26/h);  % Matriz do potencial, em X e Y

% Para cada nó, acha X e Y, e salva nas correspondentes posições da malha
len = size(nodes);
for k = 1:len(1)
  x = nodes(k,1);
  y = nodes(k,2);
  j = 1 + x/h;
  i = (y+26)/h;
  V(i,j) = Vnodes(k);
  Az(i,j) = Aznodes(k);
end

% Completa o domínio com os simétricos
V = [flip(V,2), V];
% V(1:26/h-1,:) = V(1:26/h-1,:)*10^-8;
Az = [flip(Az,2), Az];

disp("Valores em X e Y recuperados")

%% Medidas dependentes
[Ex, Ey] = gradient(-V,h); % Cálculo do campo elétrico
% Ex(1:26/h-1,:) = Ex(1:26/h-1,:)*10^-8;
% Ey(1:26/h-1,:) = Ey(1:26/h-1,:)*10^-8;

[By, Bx] = gradient(Az,h); % Cálculo do campo magnético
Bx = -Bx; % Inverte o sinal por causa do sist. de coordenadas

middle = ceil(52/(2*h));   % Posição do meio
% Pondera B por mi de cada campo, para obter H
Hx = [Bx(1:middle,:)/mi_solo;Bx(middle+1:end,:)/mi_ar];
Hy = [By(1:middle,:)/mi_solo;By(middle+1:end,:)/mi_ar];

disp("Medidas dependentes calculadas")

%% Plotagem
PlotResults(h, V, Az, Ex, Ey, Bx, By, Hx, Hy);

disp("Fim de programa")

%% Funções


function [SKV, SKB, SFV, SFB] = BoundaryConditions(h, SK, nodes, Nx, mi_ar)
  SFV = zeros(Nx,1); 
  SFB = zeros(Nx,1);
  
  SKV = SK;
  SKB = SK; 
  
  %% Fontes(cabos)
  Vmax = 500000;      % Tensão na fonte
  Jz = 200/(2*pi*0.02); % Efeito magnético na fonte
  [~, i_fonte_1] = ismember([6,10],nodes,'rows');  % Cabo de baixo
  [~, i_fonte_2] = ismember([4,14],nodes,'rows');  % Cabo de cima
  for i = 1:Nx
    SFV(i) = SFV(i) - SK(i,i_fonte_2)*Vmax - SK(i,i_fonte_1)*Vmax;
    SKV(i_fonte_1,i) = 0;
    SKV(i,i_fonte_1) = 0;
    SKV(i_fonte_2,i) = 0;
    SKV(i,i_fonte_2) = 0;
  end
  SFV(i_fonte_1) = 500000;
  SFV(i_fonte_2) = 500000;
  SKV(i_fonte_1,i_fonte_1) = 1;
  SKV(i_fonte_2,i_fonte_2) = 1;
  
  SFB(i_fonte_1) = mi_ar*Jz*(h^2/2);
  SFB(i_fonte_2) = mi_ar*Jz*(h^2/2);
  
  %% Bordas(exceto última coluna)
  for x = 0:h:26-2*h
    ylimit = sqrt(26^2 - x^2);
    innerylimit = floor((ylimit-0.00001)/h)*h;
    [~, i_topo] = ismember([x,innerylimit],nodes,'rows'); % Nó de cima, para esse x
    [~, i_baixo] = ismember([x,-innerylimit],nodes,'rows');% Nó de baixo, para esse x
    for i = 1:Nx
      SFV(i) = SFV(i) - SKV(i,i_baixo)*0 - SKV(i,i_topo)*0;
      SKV(i_topo,i) = 0;
      SKV(i,i_topo) = 0;
      SKV(i_baixo,i) = 0;
      SKV(i,i_baixo) = 0;
      
      SFB(i) = SFB(i) - SKB(i,i_baixo)*0 - SKB(i,i_topo)*0;
      SKB(i_topo,i) = 0;
      SKB(i,i_topo) = 0;
      SKB(i_baixo,i) = 0;
      SKB(i,i_baixo) = 0;
    end
    SFV(i_topo) = 0;
    SFV(i_baixo) = 0;
    SKV(i_topo,i_topo) = 1;
    SKV(i_baixo,i_baixo) = 1;
    
    SFB(i_topo) = 0;
    SFB(i_baixo) = 0;
    SKB(i_topo,i_topo) = 1;
    SKB(i_baixo,i_baixo) = 1;
  end
 
  
  %% Última coluna
  x = 26-h;
  ylimit = sqrt(26^2 - x^2);
  innerylimit = floor((ylimit-0.00001)/h)*h; 
  for y = -innerylimit:h:innerylimit % Itera todos pontos da última coluna
    [~, i_no] = ismember([x,y],nodes,'rows');
    for i = 1:Nx
      SFV(i) = SFV(i) - SKV(i,i_no)*0;
      SKV(i_no,i) = 0;
      SKV(i,i_no) = 0;
      
      SFB(i) = SFB(i) - SKB(i,i_no)*0;
      SKB(i_no,i) = 0;
      SKB(i,i_no) = 0;
    end
    SFV(i_no) = 0;
    SKV(i_no,i_no) = 1;
    
    SFB(i_no) = 0;
    SKB(i_no,i_no) = 1;
  end
  
  %% Carro
  car_x = 0:h:2;              % Lista de posições discretas em x
  car_y = 0:h:2;              % Lista de posições discretas em y
  car_y = car_y(2:end-1);     % Remove pontos extremos, para evitar repetição

  [~, width] = size(car_x);   % Obtém largura, em pontos
  [~, height] = size(car_y);  % Obtém altura, em pontos

  car = [car_x, car_x, 2*ones(1, height)];                % x de todos pontos
  car = [car; 0*ones(1,width), 2*ones(1,width), car_y];   % y de todos pontos
  car = transpose(car); % Muda orientaçao, para cada ponto ser uma linha da matriz
  
  len = size(car);
  
  % Considera ponto do lado direito como referência do carro
  x_car = car(end,1);
  y_car = car(end,2);
  [~, i_car] = ismember([x_car,y_car],nodes,'rows'); % Referência do carro
  % Itera todos pontos do carro, colocando como condição que o seu valor
  % é o mesmo que o da referência
  for k = 1:len(1)-1
    x = car(k,1);
    y = car(k,2);
    [~, i_no] = ismember([x,y],nodes,'rows');
    for i = 1:Nx
      SKV(i_no,i) = 0;
      SKB(i_no,i) = 0;
    end
    SKV(i_no,i_no) = 1;
    SKV(i_no, i_car) = -1;
    SKB(i_no,i_no) = 1;
    SKB(i_no, i_car) = -1;
  end
  
end

%% Cria Elementos
function [nodes, Nx] = CreateNodes(h)
  nodes = []; 
  Nx = 0;     
  for x = 0:h:26 
    ylimit = sqrt(26^2 - x^2);  
    innerylimit = floor((ylimit-0.00001)/h)*h;
    for y = -innerylimit:h:innerylimit 
      nodes = [nodes; x y]; 
      Nx = Nx+1;
    end
  end
end

function [K] = ElementMatrix(Ae, x1, y1, x2, y2, x3, y3)
    global sigma_ar sigma_solo mi_ar mi_solo
    % Constantes da geometria do elemento
    b1 = y2-y3;
    b2 = y3-y1;
    b3 = y1-y2;    
    c1 = x3-x2;
    c2 = x1-x3;
    c3 = x2-x1;
    
    % Matriz de "rigidez"
    K = zeros(3);
    K(1,1) = b1^2 + c1^2;
    K(1,2) = b1*b2 + c1*c2;
    K(1,3) = b1*b3 + c1*c3;  
    K(2,1) = K(1,2);
    K(2,2) = b2^2 + c2^2;
    K(2,3) = b2*b3 + c2*c3;    
    K(3,1) = K(1,3);
    K(3,2) = K(2,3);
    K(3,3) = b3^2 + c3^2;
    
    if x1 >= 0
        K = (1/(4*Ae))*K*sigma_ar;
    else 
        K = (1/(4*Ae))*K*sigma_solo;
    end
end

function [SK] = CondenseMatrix(SK, EK,i1, i2, i3)
  relationsTable = [i1 i2 i3];  
  for i = 1:3
    line = relationsTable(i);
    for j = 1:3
      column = relationsTable(j);
      SK(line,column) = SK(line, column) + EK(i,j);
    end
  end
end

function [SK] = CreateElements(h,nodes,Nx)
  SK = zeros(Nx);  
  for x = 0:h:26   
    ylimit = sqrt(26^2 - x^2); 
    innerylimit = floor((ylimit)/h)*h - h;
    
    % Metade de cima
    for y = 0:h:innerylimit 
      x1 = x;
      y1 = y;
      % Triângulo de baixo
      x2 = x+h;
      y2 = y;
      x3 = x+h;
      y3 = y+h;
      [p1, i1] = ismember([x1,y1],nodes,'rows');
      [p2, i2] = ismember([x2,y2],nodes,'rows');
      [p3, i3] = ismember([x3,y3],nodes,'rows');
      if p1 && p2 && p3
        [K] = ElementMatrix(h^2/2,x1,y1,x2,y2,x3,y3);
        [SK] = CondenseMatrix(SK,K,i1,i2,i3);
      end      
      % Triângulo de cima
      x2 = x+h;
      y2 = y+h;
      x3 = x;
      y3 = y+h;
      [p1, i1] = ismember([x1,y1],nodes,'rows');
      [p2, i2] = ismember([x2,y2],nodes,'rows');
      [p3, i3] = ismember([x3,y3],nodes,'rows');
      if p1 && p2 && p3
        [K] = ElementMatrix(h^2/2,x1,y1,x2,y2,x3,y3);
        [SK] = CondenseMatrix(SK,K,i1,i2,i3);
      end
    end
    
    % Metade de baixo
    for y = 0:-h:-innerylimit 
      x1 = x;
      y1 = y;
      % Triângulo de cima
      x2 = x+h;
      y2 = y-h;
      x3 = x+h;
      y3 = y;   
      [p1, i1] = ismember([x1,y1],nodes,'rows');
      [p2, i2] = ismember([x2,y2],nodes,'rows');
      [p3, i3] = ismember([x3,y3],nodes,'rows');   
      if p1 && p2 && p3
        [K] = ElementMatrix(h^2/2,x1,y1,x2,y2,x3,y3);
        [SK] = CondenseMatrix(SK,K,i1,i2,i3);
      end      
      % Triângulo de baixo
      x2 = x;
      y2 = y-h;
      x3 = x+h;
      y3 = y-h;
      [p1, i1] = ismember([x1,y1],nodes,'rows');
      [p2, i2] = ismember([x2,y2],nodes,'rows');
      [p3, i3] = ismember([x3,y3],nodes,'rows');
      if p1 && p2 && p3
        [K] = ElementMatrix(h^2/2,x1,y1,x2,y2,x3,y3);
        [SK] = CondenseMatrix(SK,K,i1,i2,i3);
      end
    end
  end
end

function [] = PlotResults(h, V, Az, Ex, Ey, Bx, By, Hx, Hy)
  %% Plotagem
  [X,Y] = meshgrid(-26+h:h:26, -26+h:h:26); % Malha para mostrar resultados

  % Plot das curvas de nível de V
  figure
  contour(X, Y, V, 30);
  title("V");
  axis equal

  % Plot das curvas de nível de V, com texto
  figure
  contour(X, Y, V, 'ShowText','on');
  title("V");
  axis equal

  % Plot da superfície de V
  figure
  surf(X, Y, V);
  title("V");

  % Plot das curvas de nível de Az
  figure
  contour(X, Y, Az, 30);
  title("Az");
  axis equal

  % Plot das curvas de nível de Az, com texto
  figure
  contour(X, Y, Az, 'ShowText','on');
  title("Az");
  axis equal

  % Plot da superfície de Az
  figure
  surf(X, Y, Az);
  title("Az");

  % Plot de E
  figure
  quiver(X, Y, Ex, Ey);
  title("E");
  axis equal

  % Plot de B
  figure
  quiver(X, Y, Bx, By);
  title("B");
  axis equal

  % Plot de H
  figure
  quiver(X, Y, Hx, Hy);
  title("H");
  axis equal
end

