function [nos_elementos,elementos] = criaelementos(barras)
    %Função responsável por fazer a discretização desejada, dividindo as
    %barras em elementos e armazenando tanto os elementos que surgem quanto
    %os nós.
    %Será uma matriz gigantesca dos elementos em que cada elemento será
    %salvo como
    %elementos=[x1,y1,x2,y2,di,de]
    global Elementos
    nos_elementos=[];
    elementos=[]; 

    for i=1:length(barras)
        for j=1:Elementos
            if j==1
                x1=barras(i,1);
                y1=barras(i,2);
                if length(nos_elementos)~=0
                    if length(find(nos_elementos(:,1)==x1))==0 ||...
                            length(find(nos_elementos(:,2)==y1))==0
                            nos_elementos(end+1,:)=barras(i,1:2);
                    else
                        a=(find(nos_elementos(:,1)==x1));
                        b=(find(nos_elementos(:,2)==y1));
                        jatem=0;
                        for m=1:length(a)
                            for n=1:length(b)
                                if a(m)==b(n)
                                    jatem=1;
                                end
                            end
                        end
                        if jatem==0
                            nos_elementos(end+1,:)=barras(i,1:2);
                        end
                    end
                else
                    nos_elementos(end+1,:)=barras(i,1:2);
                end
                x2=barras(i,1)+(j/Elementos)*(barras(i,3)-barras(i,1));
                y2=barras(i,2)+(j/Elementos)*(barras(i,4)-barras(i,2));
                if Elementos~=1 || i==length(barras)
                    nos_elementos(end+1,:)=[x2,y2];
                end
                elementos(end+1,:)=[barras(i,1:2),x2,y2,barras(i,5:6)];
            elseif (j==Elementos && j~=1)
                x1=x2;
                y1=y2;
                x2=barras(i,3);
                y2=barras(i,4);
                if length(find(nos_elementos(:,1)==x2))==0 ||...
                        length(find(nos_elementos(:,2)==y2))==0
                        nos_elementos(end+1,:)=barras(i,3:4);
                else
                    a=(find(nos_elementos(:,1)==x2));
                    b=(find(nos_elementos(:,2)==y2));
                    jatem=0;
                    for m=1:length(a)
                        for n=1:length(b)
                            if a(m)==b(n)
                                jatem=1;
                            end
                        end
                    end
                    if jatem==0
                        nos_elementos(end+1,:)=barras(i,3:4);
                    end
                end
                elementos(end+1,:)=[x1,y1,barras(i,3:6)];
            else
                x1=x2;
                y1=y2;
                x2=barras(i,1)+(j/Elementos)*(barras(i,3)-barras(i,1));
                y2=barras(i,2)+(j/Elementos)*(barras(i,4)-barras(i,2));
                nos_elementos(end+1,:)=[x2,y2];
                elementos(end+1,:)=[x1 y1 x2 y2 barras(i,5:6)];
            end
        end
    end

    plotatorre(nos_elementos,elementos,2,0);
end

