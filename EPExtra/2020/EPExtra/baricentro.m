function [XMedio,YMedio]=baricentro(classeelementos)
    % Função responsável por encontrar o baricentro dos pontos para aplicar
    % os vetores para evitar dois vetores aplicados no mesmo ponto
    XMedio =zeros(length(classeelementos),1);
    YMedio =zeros(length(classeelementos),1);
    for i = 1:length(classeelementos)
        XMedio(i)=(classeelementos(i).no1posx+classeelementos(i).no2posx+...
            classeelementos(i).no3posx)/3;
        YMedio(i)=(classeelementos(i).no1posy+classeelementos(i).no2posy+...
            classeelementos(i).no3posy)/3;    
    end
end

