function Item1(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g,Theta1_Inicial,Theta2_Inicial,Alpha1_Inicial,Alpha2_Inicial,h,t_Inicial,t_Final)
%Item1 = Realiza o cálculo da dinâmica do sistema utilizando os 3 métodos
%para 6 diferentes valores de h, analisando assim a influência tanto do
%método quanto do passo.

%Gera o vetor de funções anônimas que são as derivadas
Derivadas=geraDerivadas(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g);

for i=1:length(h);
    tempo=t_Inicial:h(i):t_Final;

    Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial];

    %Roda e plota o método de Euler
    [Variaveis, Arquiva_Derivadas]=Euler(Variaveis,Derivadas,h(i),tempo);
    defineGraficos(tempo,Variaveis,Arquiva_Derivadas,1,h(i),1,1);

    %Roda e plota o método RK2
    Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial]; 
    %Redefine o vetor Variaveis para não influenciar o outro método

    [Variaveis, Arquiva_Derivadas]=RK2(Variaveis,Derivadas,h(i),tempo); 
    defineGraficos(tempo,Variaveis,Arquiva_Derivadas,2,h(i),1,1);

    %Roda e plota o método RK4
    Variaveis=[Theta1_Inicial;Theta2_Inicial;Alpha1_Inicial;Alpha2_Inicial]; 

    [Variaveis, Arquiva_Derivadas]=RK4(Variaveis,Derivadas,h(i),tempo);
    defineGraficos(tempo,Variaveis,Arquiva_Derivadas,3,h(i),1,1);
end
end

