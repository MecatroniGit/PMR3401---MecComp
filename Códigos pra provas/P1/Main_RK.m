%Variáveis iniciais
%Aqui, define tudo que é relevante como variável inicial
%Importante lembrar que as variáveis são as variáveis DEPENDENTES, a
%variável em relação a qual elas são derivadas NÃO ENTRA AQUI
%Variaveis=[y0,z0,w0]
VariaveisIniciais=[0;2];
Variaveis=VariaveisIniciais;

%Aqui define quantos passos você vai dar e instantes
h=0.1;
Tinicial=0;
Tfinal=0.2;
%Chamamos de TEMPO inicial e final, mas caso esteja sendo derivado em
%relação a X, será X inicial e final, não precisa mudar as variáveis.
tempo=Tinicial:h:Tfinal;

%Define qual RungeKutta você quer
% Euler - 1
% RK2 - Heun - 1
% RK2 - Euler Modificado - 2
% RK4 - 4
RK=2;

%Deseja calcular o erro de passo dobrado?
%Sim - 1
%Não - 0
Desejo=1;

%Define os valores pra por no vetor das derivadas (que vão ser usados pra
%calcular K1
Derivadas=derivadas();
if Desejo==1
    VariaveisFinais1=[0 0]
    for i=1:2
        if RK==0
            Variaveis=Euler(Variaveis,Derivadas,h,tempo);
        elseif RK~=4
            Variaveis=RK2(Variaveis,Derivadas,h,tempo,RK);
        else
            Variaveis=RK4(Variaveis,Derivadas,h,tempo);
        end
        VariaveisFinais(:,i)=Variaveis(:,end)
        h=h/2
        tempo=Tinicial:h:Tfinal;
        Variaveis=VariaveisIniciais;
    end
    if RK~4 && RK~=0
        m=2;
    else
        m=4;
    end
    erro=(2^m*(VariaveisFinais(:,2)-VariaveisFinais(:,1)))/(2^m-1)
else
    if RK==0
        Variaveis=Euler(Variaveis,Derivadas,h,tempo);
    elseif RK~=4
        Variaveis=RK2(Variaveis,Derivadas,h,tempo,RK);
    else
        Variaveis=RK4(Variaveis,Derivadas,h,tempo);
    end
end
