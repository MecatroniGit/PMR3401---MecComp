%Vari�veis iniciais
%Aqui, define tudo que � relevante como vari�vel inicial
%Importante lembrar que as vari�veis s�o as vari�veis DEPENDENTES, a
%vari�vel em rela��o a qual elas s�o derivadas N�O ENTRA AQUI
%Variaveis=[y0,z0,w0]
VariaveisIniciais=[0;2];
Variaveis=VariaveisIniciais;

%Aqui define quantos passos voc� vai dar e instantes
h=0.1;
Tinicial=0;
Tfinal=0.2;
%Chamamos de TEMPO inicial e final, mas caso esteja sendo derivado em
%rela��o a X, ser� X inicial e final, n�o precisa mudar as vari�veis.
tempo=Tinicial:h:Tfinal;

%Define qual RungeKutta voc� quer
% Euler - 1
% RK2 - Heun - 1
% RK2 - Euler Modificado - 2
% RK4 - 4
RK=2;

%Deseja calcular o erro de passo dobrado?
%Sim - 1
%N�o - 0
Desejo=1;

%Define os valores pra por no vetor das derivadas (que v�o ser usados pra
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
