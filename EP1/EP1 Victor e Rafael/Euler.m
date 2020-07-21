function [Variaveis, Arquiva_Derivadas]= Euler(Variaveis,Derivadas,h,tempo)
%Euler = Fun��o para executar o m�todo de Euler para um sistema de EDOs
%gen�rico.

for n=1:length(tempo)-1;
    K1=Derivadas(tempo(n),Variaveis(:,n))
    Variaveis(:,n+1)= Variaveis(:,n)+h*K1
end

end
