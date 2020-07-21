function [Klocal] = montalocal(classeelementos)
    % Funcao responsavel por montar o K de cada elemento
    [Klocal,bc]=prodinterno(classeelementos);
    area=det(bc(1:2,1:2))/2;
    Klocal=Klocal/(4*area);
end