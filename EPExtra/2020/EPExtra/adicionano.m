function [listanos,classeelementos]=...
    adicionano(no,listanos,classeelementos,i,qualno)
    % Fun��o respons�vel por criar a lista de n�s e atualizar os elementos
    % com a posi��o dos n�s aos quais est�o presos
    if qualno==1
        if isempty(no)
            listanos(end+1,:)=[classeelementos(i).no1posx...
                classeelementos(i).no1posy];
            classeelementos(i).no1elemento=length(listanos(:,1));
        else
            classeelementos(i).no1elemento=no(1);
        end
    elseif qualno==2
        if isempty(no)
            listanos(end+1,:)=[classeelementos(i).no2posx...
                classeelementos(i).no2posy];
            classeelementos(i).no2elemento=length(listanos(:,1));
        else
            classeelementos(i).no2elemento=no(1);
        end
    else
        if isempty(no)
            listanos(end+1,:)=[classeelementos(i).no3posx...
                classeelementos(i).no3posy];
            classeelementos(i).no3elemento=length(listanos(:,1));
        else
            classeelementos(i).no3elemento=no(1);
        end    
    end
end