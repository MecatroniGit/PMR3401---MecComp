function [Variaveis, Arquiva_Derivadas]= RK4(Variaveis,Derivadas,h,tempo)
%RK4 = Função para executar o método RK4 para um sistema de EDOs
%genérico.

Arquiva_Derivadas=[0.4;-0.1;0;0];
for n=1:length(tempo)-1;
    K1=Derivadas(tempo(n),Variaveis(:,n));
    K2=Derivadas(tempo(n)+h/2,Variaveis(:,n)+h/2*K1);
    K3=Derivadas(tempo(n)+h/2,Variaveis(:,n)+h/2*K2);
    K4=Derivadas(tempo(n)+h,Variaveis(:,n)+h*K3);
    K=K1+2*K2+2*K3+K4;
    Arquiva_Derivadas(:,n+1)=K1;
    Variaveis(:,n+1)= Variaveis(:,n)+h/6*K;
end

end

