function Variaveis=RK4(Variaveis,Derivadas,h,tempo)
%RK4 = Função para executar o método RK4 para um sistema de EDOs
%genérico.

for n=1:length(tempo)-1;
    display('O passo em que estamos é:')
    passo=tempo(n)/h
    display('E os valores dos Ks são:')
    display(['K1=F(',num2str(tempo(n)),' ',num2str(Variaveis(:,n).'), ')'])
    K1=Derivadas(tempo(n),Variaveis(:,n))
    display(['K2=F(',num2str(tempo(n)+0.5*h),' ',num2str((Variaveis(:,n)+0.5*h*K1).'), ')'])
    K2=Derivadas(tempo(n)+h/2,Variaveis(:,n)+h/2*K1)
    display(['K3=F(',num2str(tempo(n)+0.5*h),' ',num2str((Variaveis(:,n)+0.5*h*K2).'), ')'])
    K3=Derivadas(tempo(n)+h/2,Variaveis(:,n)+h/2*K2)
    display(['K4=F(',num2str(tempo(n)+h),' ',num2str((Variaveis(:,n)+h*K3).'), ')'])
    K4=Derivadas(tempo(n)+h,Variaveis(:,n)+h*K3)
    K=K1+2*K2+2*K3+K4
    display('E as variáveis ao fim são:')
    Variaveis(:,n+1)= Variaveis(:,n)+h/6*K;
    Variaveizinhas=Variaveis(:,length(Variaveis(1,:)))
end

end

