function [pressaoarrumada,pontosoriginais]=ItemD(psi,H,h,d,L,dx,dy,pressaovetor,bordairregular,AB)
%Resolve o Item D da Parte 1 do EP2

pressao=[0 0 0];
pontosoriginais=[0 0];
for i=2:length(bordairregular)
    pontosoriginais(end+1,:)=[dx*(bordairregular(i,2)-1),H-dy*(bordairregular(i,1)-1)];
    if AB(i,1)>1
        if AB(i,2)<=1
            pressao(end+1,:)=[dx*(bordairregular(i,2)-1),H-dy*(bordairregular(i,1)-1)-AB(i,2)*dy,pressaovetor(bordairregular(i,1),bordairregular(i,2))];
        else
            %Só vai entrar aqui se for o ponto do topo que defini com
            %distâncias (2,2), e esse não precisa subtrair nada
            pressao(end+1,:)=[dx*(bordairregular(i,2)-1),H-dy*...
                (bordairregular(i,1)-1),pressaovetor(bordairregular(i,1)...
                ,bordairregular(i,2))];
        end
    elseif AB(i,2)>1
        if AB(i,1)<=1
            if bordairregular(i,2)>(length(psi(1,:))-1)/2
                pressao(end+1,:)=[dx*(bordairregular(i,2)-1)-AB(i,1)*...
                    dx,H-dy*(bordairregular(i,1)-1),...
                    pressaovetor(bordairregular(i,1),bordairregular(i,2))];
            else
                pressao(end+1,:)=[dx*(bordairregular(i,2)-1)+AB(i,1)*dx,...
                    H-dy*(bordairregular(i,1)-1),...
                    pressaovetor(bordairregular(i,1),bordairregular(i,2))];    
            end
        end
    else
        if AB(i,1)>=AB(i,2)
            %Ajusta a distancia de acordo com o valor de b, caso seja menor
            %que a
            pressao(end+1,:)=[dx*(bordairregular(i,2)-1),H-dy*...
                (bordairregular(i,1)-1)-AB(i,2)*dy,...
                pressaovetor(bordairregular(i,1),bordairregular(i,2))];
        else
            %Ajusta a distancia de acordo com o valor de a, caso seja menor
            %que b
            if bordairregular(i,2)>(length(psi(1,:))-1)/2
                pressao(end+1,:)=[dx*(bordairregular(i,2)-1)-AB(i,1)*dx,...
                    H-dy*(bordairregular(i,1)-1),...
                    pressaovetor(bordairregular(i,1),bordairregular(i,2))];
            else
                pressao(end+1,:)=[dx*(bordairregular(i,2)-1)+AB(i,1)*dx,...
                    H-dy*(bordairregular(i,1)-1),...
                    pressaovetor(bordairregular(i,1),bordairregular(i,2))];    
            end                
        end
    end
end

pressao=[pressao pontosoriginais];
[arrumado,indices]=sort(pressao(:,1));
pressaoarrumada=pressao(indices,1:3);
pontosoriginais=pressao(indices,4:5);

end

