function [pontosborda,pontosfonte,pontoscarro]=geraCCs(listanos)
    % Função responsável por encontrar os pontos nos quais temos que
    % aplicar as condições de contorno
    global raio h borda fonte1 fonte2
    pontosborda=find(sqrt(listanos(:,1).^2+listanos(:,2).^2)>raio-h-borda);
    pontosfonte=find(listanos(:,1)==fonte1(1) & listanos(:,2)==fonte1(2));
    pontosfonte(end+1)=find(listanos(:,1)==fonte2(1) & ...
        listanos(:,2)==fonte2(2));
    pontoscarro=find(listanos(:,1)<=2 &...
        listanos(:,2)<=1.5 & listanos(:,2)>=0);
end

