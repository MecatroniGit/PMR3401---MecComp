function classeelementos=criaelementos
    % Função responsável por criar a classe elementos
    global h raio
    classeelementos = [];
    for  y=0:h:raio-h
        x=0;
        limite=sqrt((raio-h)^2 - y^2);    
        while x<= limite
            x1=x;
            x2=x+h;
            y1=y;
            y2=y+h;       
            classeelementos=[classeelementos...
                montaelemento(x1,x2,x2,y1,y1,y2)];
            classeelementos=[classeelementos...
                montaelemento(x1,x2,x1,y1,y2,y2)];
            classeelementos=[classeelementos...
                montaelemento(x1,x2,x2,-y1,-y2,-y1)];
            classeelementos=[classeelementos...
                montaelemento(x1,x1,x2,-y1,-y2,-y2)];
            x = x+h; 
        end
    end
end

