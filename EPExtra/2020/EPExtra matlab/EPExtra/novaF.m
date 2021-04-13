function [Forca,Kglobal]=novaF(listanos,pontosborda,pontosfonte,...
    pontoscarro,Kglobal,parametro)
    % Função responsável por gerar o vetor de carregamentos Forca
    global mi1 Imaximo RC Vmaximo

    Kglobal(pontosborda,:)=0;
    Kglobal(:,pontosborda)=0;
    Kglobal(pontosborda,pontosborda)=Kglobal(pontosborda,pontosborda)+...
        diag(1+diag(Kglobal(pontosborda,pontosborda)));
    
    if parametro==1
        Kglobal(pontoscarro,:)=0;
        Kglobal(:,pontoscarro)=0;
        Kglobal(pontoscarro,pontoscarro)=Kglobal(pontoscarro,pontoscarro)+...
            diag(1+diag(Kglobal(pontoscarro,pontoscarro)));
    end

    Forca = zeros(length(listanos),1);

    for i = 1:length(pontosfonte)
        if parametro==1
            Forca=Forca-Kglobal(:,pontosfonte(i))*Vmaximo;
            Forca(pontosfonte(i))=Vmaximo;
        else
            Forca=Forca-Kglobal(:,pontosfonte(i))*-mi1*Imaximo/(2*pi*RC);
            Forca(pontosfonte(i))=-mi1*Imaximo/(2*pi*RC);
        end
        Kglobal(pontosfonte(i),:)=0;
        Kglobal(:,pontosfonte(i))=0;
        Kglobal(pontosfonte(i),pontosfonte(i))=1;
    end

end

