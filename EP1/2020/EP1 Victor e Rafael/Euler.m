function [Variaveis, Arquiva_Derivadas]= Euler(Variaveis,Derivadas,h,tempo)
%Euler = Função para executar o método de Euler para um sistema de EDOs
%genérico.

for n=1:length(tempo)-1;
    K1=Derivadas(tempo(n),Variaveis(:,n))
    Variaveis(:,n+1)= Variaveis(:,n)+h*K1
end

end
