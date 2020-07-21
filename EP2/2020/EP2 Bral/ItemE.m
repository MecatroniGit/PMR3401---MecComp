function [FR,Areas,Normais,Forca]=ItemE(pressaoarrumada,Comprimento)
%Resolve o Item E da Parte 1 do EP2

Normais=zeros(([length(pressaoarrumada)-1 2]));
Pontos_medios=zeros(([length(pressaoarrumada) 2]));

for i=2:(length(pressaoarrumada))
    if i==2 
        Normais(i-1,:)=[3-pressaoarrumada(i+1,2) pressaoarrumada(i+1,1)-...
            15];
        Pontos_medios(i-1,:)=[(pressaoarrumada(i,1)+15)/2 (3+...
            pressaoarrumada(i,2))/2];
    elseif i==length(pressaoarrumada)
        Normais(i-1,:)=[pressaoarrumada(i-1,2)-3 21-...
            pressaoarrumada(i-1,1)];
        Pontos_medios(i-1,:)=[(pressaoarrumada(i-1,1)+...
        pressaoarrumada(i,1))/2 (pressaoarrumada(i,2)+...
        pressaoarrumada(i-1,2))/2 ];
    else
        %Gera os vetores que ligam os pontos das press�es em i+1 e i-1,
        %que teria coordenadas (x,y) e salva em normais como (-y,x), o que
        %� um vetor normal a esse anterior.
        Normais(i-1,:)=[pressaoarrumada(i-1,2)-pressaoarrumada(i+1,2),...
            pressaoarrumada(i+1,1)-pressaoarrumada(i-1,1)];
        Pontos_medios(i-1,:)=[(pressaoarrumada(i-1,1)+...
        pressaoarrumada(i,1))/2 (pressaoarrumada(i-1,2)+...
        pressaoarrumada(i,2))/2];
    end
    Normais(i-1,:)=Normais(i-1,:)/sqrt(Normais(i-1,1)^2+Normais(i-1,2)^2);
end

Pontos_medios(end,:)=[(pressaoarrumada(end,1)+21)/2 (3+...
    pressaoarrumada(end,2))/2];
%Calcula as �reas necess�rias no EP
Areas=Comprimento*(sqrt((Pontos_medios(2:end,1)-...
    Pontos_medios(1:end-1,1)).^2+(Pontos_medios(2:end,2)-...
    Pontos_medios(1:end-1,2)).^2));
Areas(1)=Areas(1)+Comprimento*(sqrt((15-Pontos_medios(1,1)).^2+...
    (3-Pontos_medios(1,2)).^2));
Areas(end)=Areas(end)+Comprimento*(sqrt((15-Pontos_medios(1,1)).^2+...
    (3-Pontos_medios(1,2)).^2));
%Calcula a for�a
Forca=Areas.*pressaoarrumada(2:end,3);
Componentey=-Forca.*Normais(:,2);
%Soma as componentes da for�a em y para obter a for�a resultante
FR=sum(Componentey);

end

