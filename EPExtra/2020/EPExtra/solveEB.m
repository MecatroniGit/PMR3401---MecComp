function [Vetx,Vety] = solveEB(classeelementos, vetor,item)
    % Função responsável por encontrar os valores de E e B como pós
    % processamento após encontrar V e Az, respectivamente.
    Vetx=zeros(length(classeelementos),1);
    Vety=zeros(length(classeelementos),1);
    for i = 1:length(classeelementos)
        element=classeelementos(i);
        if item==1
            vet=[vetor(element.no1elemento) vetor(element.no2elemento)...
                vetor(element.no3elemento)];
        else
            vet=[vetor(element.no1elemento) vetor(element.no2elemento)...
                vetor(element.no3elemento)];
        end
        [Klocal,bc]=prodinterno(element);
        area=det(bc(1:2,1:2))/2;
        dx=(1/(2*area))*dot(bc(:,1)',vet);
        dy=(1/(2*area))*dot(bc(:,2)',vet);
        if item==1
            Vetx(i) = -dx;
            Vety(i) = -dy;
        elseif item==2
            Vetx(i)=dy;
            Vety(i)=-dx;
        end
    end
      
end
