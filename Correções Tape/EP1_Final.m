%Constantes do problema
L1=2;
L2=2.5;
L2eixo=1.8;
g=9.81;
m1=450;
m2=650;
F1=-0.5*m1*g;
F2=-0.5*m2*g;
miIz=2.7;
R=0.3;
Vel=80/3.6;

%Tempo e passo
t_Inicial=0;
t_Final=60;
h= [0.1 0.05 0.01 0.005 0.001];

%Condições iniciais
Alpha1_Inicial=0.4; %Theta ponto 1
Theta1_Inicial=0; 
Alpha2_Inicial=-0.1; %Theta ponto 2
Theta2_Inicial=0;

%Item 1: Comparar métodos e influência do passo h nas respostas
Item1(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g,Theta1_Inicial,Theta2_Inicial,Alpha1_Inicial,Alpha2_Inicial,h,t_Inicial,t_Final)

%Item 2: Comparar influência de diferentes condições iniciais na dinâmica
%do sistema
Item2(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g,Theta1_Inicial,Theta2_Inicial,Alpha1_Inicial,Alpha2_Inicial,h,t_Inicial,t_Final)