function barras = criabarras(nos)
    %Função responsável por criar as matrizes que armazenam as barras da
    %torre, posteriormente discretizadas em elementos
    %A organização dentro da matriz é de cada linha ser uma barra e os
    %valores contidos dentro serem
    %barras=[x1 y1 x2 y2 di de]
    global d1i d1e d2i d2e
    barras=zeros(27,6);
    y=1;
    direita=1;
    for i=1:2:length(nos)-2
        if nos(i,1)<nos(i+2,1)
            if nos(i,1)==0
                barras(y,:)=...
                    [nos(i,1),nos(i,2),nos(i+1,1),nos(i+1,2),d2i,d2e];
            else
                barras(y,:)=...
                    [nos(i,1),nos(i,2),nos(i+1,1),nos(i+1,2),d1i,d1e];                
            end
            barras(y+1,:)=...
                [nos(i,1),nos(i,2),nos(i+2,1),nos(i+2,2),d2i,d2e];
            barras(y+2,:)=...
                [nos(i+1,1),nos(i+1,2),nos(i+3,1),nos(i+3,2),d2i,d2e];
            barras(y+3,:)=...
                [nos(i,1),nos(i,2),nos(i+3,1),nos(i+3,2),d1i,d1e];
            barras(y+4,:)=...
                [nos(i+1,1),nos(i+1,2),nos(i+2,1),nos(i+2,2),d1i,d1e];
            y=y+5;
        else
            barras(y,:)=[nos(i,1),nos(i,2),nos(i+1,1),nos(i+1,2),d1i,d1e];
            barras(y+1,:)=...
                [nos(i,1),nos(i,2),nos(i+2,1),nos(i+2,2),d2i,d2e];
            barras(y+2,:)=...
                [nos(i+1,1),nos(i+1,2),nos(i+3,1),nos(i+3,2),d2i,d2e];
            if direita==1
                barras(y+3,:)=...
                    [nos(i,1),nos(i,2),nos(i+3,1),nos(i+3,2),d1i,d1e];
                direita=0;
            else
                barras(y+3,:)=...
                    [nos(i+1,1),nos(i+1,2),nos(i+2,1),nos(i+2,2),d1i,d1e];
                direita=1;
            end
            y=y+4;
        end
    end
    barras(y,:)=[nos(length(nos)-1,1),...
      nos(length(nos)-1,2),nos(length(nos),1),nos(length(nos),2),d2i,d2e];
    plotatorre(nos,barras,1,0);
end

