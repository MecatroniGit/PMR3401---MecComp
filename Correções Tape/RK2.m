function RK2(Variaveis,Derivadas,h,tempo,RK)
%RK2 = Fun��o para executar o m�todo RK2 para um sistema de EDOs
%gen�rico.

if RK==1
    for n=1:length(tempo)-1;
        printf('O instante de tempo em que estamos �:')
        tempo(n)
        printf('E os valores dos Ks s�o:')
        K1=Derivadas(tempo(n),Variaveis(:,n))
        K2=Derivadas(tempo(n)+h,Variaveis(:,n)+h*K1)
        K=K1+K2;
        printf('E as vari�veis ao fim s�o:')
        Variaveis(:,n+1)= Variaveis(:,n)+h/2*K;
    end
else
    for n=1:length(tempo)-1;
        printf('O instante de tempo em que estamos �:')
        tempo(n)
        printf('E os valores dos Ks s�o:')
        K1=Derivadas(tempo(n),Variaveis(:,n))
        K2=Derivadas(tempo(n)+0.5*h,Variaveis(:,n)+0.5*h*K1)
        printf('E as vari�veis ao fim s�o:')
        Variaveis(:,n+1)= Variaveis(:,n)+h*K2;
    end    
end
end

