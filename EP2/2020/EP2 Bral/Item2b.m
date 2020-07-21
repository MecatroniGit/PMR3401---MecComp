function [calor,gradientesarrumado,Normais]=Item2b(T,bordairregular,AB,psi,dx,dy,H,h,L,d,Areas,Normais,k)
%Resolve o Item B da Parte 2 do EP2

gradientes=[0 0 0 0]; %x, y, dT/dx dT/dy
for i=2:length(bordairregular)
    if bordairregular(i,2)>(length(psi(1,:))+1)/2
        %u>0 e v>0
        if AB(i,1)<1
            dTdx=(T(bordairregular(i,1),bordairregular(i,2)-1)-...
            T(bordairregular(i,1),bordairregular(i,2)))/dx;
        else
            dTdx=(T(bordairregular(i,1),bordairregular(i,2)-1)-...
            T(bordairregular(i,1),bordairregular(i,2)))/dx;          
        end
        dTdy=(T(bordairregular(i,1),bordairregular(i,2))-...
        T(bordairregular(i,1)-1,bordairregular(i,2)))/dx;
        if AB(i,1)>=AB(i,2)
            gradientes(end+1,:)=[dx*(bordairregular(i,2)-1),...
            H-dy*(bordairregular(i,1)-1)-AB(i,2)*dy,dTdx,dTdy];
        else
            gradientes(end+1,:)=[dx*(bordairregular(i,2)-1)-AB(i,1)*dx,...
            H-dy*(bordairregular(i,1)-1),dTdx,dTdy];
        end
    elseif bordairregular(i,2)<(length(psi(1,:))+1)/2
        %u>0 e v<0
        if AB(i,1)<1
            dTdx=(T(bordairregular(i,1),bordairregular(i,2))-...
            T(bordairregular(i,1),bordairregular(i,2)-1))/dx;
        else
            dTdx=(T(bordairregular(i,1),bordairregular(i,2))-...
            T(bordairregular(i,1),bordairregular(i,2)-1))/dx;            
        end
        if AB(i,2)<1
            dTdy=(T(bordairregular(i,1)+1,bordairregular(i,2))-...
            T(bordairregular(i,1),bordairregular(i,2)))/dx;
        else
            dTdy=(T(bordairregular(i,1)+1,bordairregular(i,2))-...
            T(bordairregular(i,1),bordairregular(i,2)))/dx;    
        end
        if AB(i,1)>=AB(i,2)
            gradientes(end+1,:)=[dx*(bordairregular(i,2)-1),...
            H-dy*(bordairregular(i,1)-1)-AB(i,2)*dy,dTdx,dTdy];
        else
            gradientes(end+1,:)=[dx*(bordairregular(i,2)-1)+AB(i,1)*dx,...
            H-dy*(bordairregular(i,1)-1),dTdx,dTdy];
        end        
    else
        %Está bem no topo
        dTdy=(T(bordairregular(i,1)-2,bordairregular(i,2))-...
        4*T(bordairregular(i,1)-1,bordairregular(i,2))+...
        3*T(bordairregular(i,2),bordairregular(i,1)))/(2*dx);   
        dTdx=(T(bordairregular(i,1),bordairregular(i,2))-...
        T(bordairregular(i,1),bordairregular(i,2)-1))/dx;            
        gradientes(end+1,:)=[dx*(bordairregular(i,2)-1),...
        H-dy*(bordairregular(i,1)-1),dTdx,dTdy];
    end
end

for j=(length(psi(:,1))-1)+1-3/dy:(length(psi(:,1)))
    dTdx=(T(j,15/dx+1)-T(j,15/dx))/dx;
    dTdy=0;
    gradientes(end+1,:)=[15,H-dx*(j-1),dTdx,dTdy];
    dTdx=(T(j,21/dx+1)-T(j,21/dx+2))/dx;
    Normais=[[-1 0];Normais];
    gradientes(end+1,:)=[21,H-dx*(j-1),dTdx,dTdy];
    Normais(end+1,:)=[1 0];
end

[arrumado,indices]=sort(gradientes(:,1));
gradientesarrumado=gradientes(indices,:);


produto=gradientesarrumado(2:end,3:4).*Normais(:,1:2);
modulo=sqrt(produto(:,1).^2+produto(:,2).^2);

lateral=zeros(3/dx-1,1);
lateral(:,1)=dx*60; 
%acrescenta nas áreas a parcela relativa aos lados do silo
Areas=[0.5*dx*60;lateral;0.5*dx*60;Areas;0.5*dx*60;lateral;0.5*dx*60];
calor=-k*sum(Areas.*modulo);

end

