function Variaveis=Euler(Variaveis,Derivadas,h,tempo)
%Euler = Função para executar o método de Euler para um sistema de EDOs
%genérico.

for n=1:length(tempo)-1;
    display('O passo em que estamos é:')
    passo=tempo(n)/n
    display(['F(',num2str(tempo(n)),' ',num2str(Variaveis(:,n).'), ')'])
    K1=Derivadas(tempo(n),Variaveis(:,n))
    display('E as variáveis finais são:')
    Variaveis(:,n+1)= Variaveis(:,n)+h*K1;
    Variaveizinhas=Variaveis(:,length(Variaveis(1,:)))
end

end
