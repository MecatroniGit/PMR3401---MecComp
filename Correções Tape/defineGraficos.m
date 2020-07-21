function defineGraficos(tempo,Variaveis,Arquiva_Derivadas,metodo,passo,item,subitem)
%defineGraficos = Identifica qual � o gr�fico que queremos plotar, define o
%t�tulo dele e chama a fun��o que realiza o plot.

passo=num2str(passo);

if item==1;
    %Define o m�todo que estamos rodando no momento
    if metodo==1;    
        titulo=strcat("Resultado com m�todo de Euler com passo ",passo);
    elseif metodo==2;
        titulo=strcat("Resultado com RK2 com passo ",passo);
    elseif metodo==3;
        titulo=strcat("Resultado com RK4 com passo ",passo);
    end

elseif item==2;
    %Define o m�todo que estamos rodando no momento
    if subitem==1;    
        titulo=strcat("Resultado com m2=1000kg e passo ",passo);
    elseif subitem==2;
        titulo=strcat("Resultado com m2=200kg e passo ",passo);        
    elseif subitem==3;
        titulo=strcat("Resultado com Vel=120km/h e passo ",passo);        
    elseif subitem==4;
        titulo=strcat("Resultado com F1=0,5*m1*g N e passo ",passo);        
    end
    
    
end

plotGraficos(titulo,tempo,Variaveis,Arquiva_Derivadas);

end
