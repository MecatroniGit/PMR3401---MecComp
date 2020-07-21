function [listanos,classeelementos]=crianos(classeelementos)
    % Função responsável por criar a lista de nós
    listanos=[];
    for i = 1:length(classeelementos)
        if isempty(listanos)
            no1=[];
        else
            no1=find(listanos(:,1)==classeelementos(i).no1posx &...
                listanos(:,2)==classeelementos(i).no1posy);
        end
        [listanos,classeelementos]=...
            adicionano(no1,listanos,classeelementos,i,1);
        no2=find(listanos(:,1)==classeelementos(i).no2posx &...
            listanos(:,2)==classeelementos(i).no2posy);    
        no3=find(listanos(:,1)==classeelementos(i).no3posx &...
            listanos(:,2)==classeelementos(i).no3posy);
        [listanos,classeelementos]=...
            adicionano(no2,listanos,classeelementos,i,2);
        [listanos,classeelementos]=...
            adicionano(no3,listanos,classeelementos,i,3);
    end
    plotadiscretizacao(listanos,classeelementos)

end

