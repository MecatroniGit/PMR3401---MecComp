function Item2(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g,Theta1_Inicial,Theta2_Inicial,Alpha1_Inicial,Alpha2_Inicial,h,t_Inicial,t_Final)
%Item 2 = Compara como a dinâmica do sistema varia com a modificação das
%condições iniciais.

tempo=t_Inicial:h(end):t_Final;

%Item 2a - aumento na massa
Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial];
m2=1000;
Derivadas=geraDerivadas(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g);

[Variaveis, Arquiva_Derivadas]=RK4(Variaveis,Derivadas,h(end),tempo);
defineGraficos(tempo,Variaveis,Arquiva_Derivadas,3,h(end),2,1);

%Item 2b - diminuição na massa
Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial];
m2=200;
Derivadas=geraDerivadas(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g);

[Variaveis, Arquiva_Derivadas]=RK4(Variaveis,Derivadas,h(end),tempo);
defineGraficos(tempo,Variaveis,Arquiva_Derivadas,3,h(end),2,2);

%Item 2c - Aumento na velocidade
Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial];
m2=650;
Vel=120/3.6;
Derivadas=geraDerivadas(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g);

[Variaveis, Arquiva_Derivadas]=RK4(Variaveis,Derivadas,h(end),tempo);
defineGraficos(tempo,Variaveis,Arquiva_Derivadas,3,h(end),2,3);

%Item 2d - Tração traseira
Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial];
Vel=80/3.6;
F1=0.5*m1*g;
Derivadas=geraDerivadas(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g);

[Variaveis, Arquiva_Derivadas]=RK4(Variaveis,Derivadas,h(end),tempo);
defineGraficos(tempo,Variaveis,Arquiva_Derivadas,3,h(end),2,4);

end

