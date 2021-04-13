function K=criaglobal(listanos,classeelementos,parametro)
    % Função responsável por criar a matriz global K para os casos de V e
    % de Az
    global mi1 mi2 sigma1 sigma2
    if parametro==1
        parametro1=sigma1;
        parametro2=sigma2;
    else
        parametro1=mi1;
        parametro2=mi2;
    end
    K = zeros(length(listanos));
    for i = 1:length(classeelementos)    
        Klocal=montalocal(classeelementos(i));
        atuzalista=[classeelementos(i).no1elemento...
         classeelementos(i).no2elemento classeelementos(i).no3elemento];
        if classeelementos(i).no3posy>=0 &&...
            classeelementos(i).no2posy>=0 && classeelementos(i).no1posy>=0
            Klocal=parametro1*Klocal;
        else
            Klocal=parametro2*Klocal;
        end
        K(atuzalista,atuzalista)=K(atuzalista,atuzalista)+Klocal;  
    end
end

