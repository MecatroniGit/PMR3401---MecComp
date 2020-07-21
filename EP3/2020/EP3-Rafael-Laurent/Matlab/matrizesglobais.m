function [Kglobs,Mglobs,Mglobs2,Lele]=...
    matrizesglobais(elementos,nos_elementos,K,M,Malternas)

%Função para, uma vez criar as matrizes globais necessárias
Mglobs=zeros(3*length(nos_elementos));
Kglobs=zeros(3*length(nos_elementos));
Mglobs2=zeros(3*length(nos_elementos));

Lele=zeros(1,length(elementos));
Iele=zeros(1,length(elementos));
Aele=zeros(1,length(elementos));

for i=1:length(elementos)
    
    Lele(i)=sqrt((elementos(i,1)-elementos(i,3))^2+...
        (elementos(i,2)-elementos(i,4))^2);
    Iele(i)=pi*(elementos(i,6)^4-elementos(i,5)^4)/64;
    Aele(i)=pi*((elementos(i,6)/2)^2-(elementos(i,5)/2)^2);
    cosseno=(elementos(i,3)-elementos(i,1))/Lele(i);
    seno=(elementos(i,4)-elementos(i,2))/Lele(i);
    T=[cosseno seno 0 0 0 0; -seno cosseno 0 0 0 0; 0 0 1 0 0 0;...
        0 0 0 cosseno seno 0; 0 0 0 -seno cosseno 0; 0 0 0 0 0 1];
    Mele=transpose(T)*M(Lele(i),Aele(i))*(T);
    
    Mele2=transpose(T)*Malternas(Lele(i),Aele(i))*(T);
    
    Kele=transpose(T)*K(Lele(i),Aele(i),Iele(i))*(T);
    achax1=find(nos_elementos(:,1)==elementos(i,1));
    achay1=find(nos_elementos(achax1,2)==elementos(i,2));
    achax2=find(nos_elementos(:,1)==elementos(i,3));
    achay2=find(nos_elementos(achax2,2)==elementos(i,4));
    a=3*achax1(achay1);
    b=3*achax2(achay2);
    Mglobs(a-2:a,a-2:a)=Mglobs(a-2:a,a-2:a)+Mele(1:3,1:3);
    Mglobs(b-2:b,b-2:b)=Mglobs(b-2:b,b-2:b)+Mele(4:6,4:6);
    Mglobs(a-2:a,b-2:b)=Mglobs(a-2:a,b-2:b)+Mele(1:3,4:6);
    Mglobs(b-2:b,a-2:a)=Mglobs(b-2:b,a-2:a)+Mele(4:6,1:3);

    Mglobs2(a-2:a,a-2:a)=Mglobs2(a-2:a,a-2:a)+Mele2(1:3,1:3);
    Mglobs2(b-2:b,b-2:b)=Mglobs2(b-2:b,b-2:b)+Mele2(4:6,4:6);
    Mglobs2(a-2:a,b-2:b)=Mglobs2(a-2:a,b-2:b)+Mele2(1:3,4:6);
    Mglobs2(b-2:b,a-2:a)=Mglobs2(b-2:b,a-2:a)+Mele2(4:6,1:3);    
    
    Kglobs(a-2:a,a-2:a)=Kglobs(a-2:a,a-2:a)+Kele(1:3,1:3);
    Kglobs(b-2:b,b-2:b)=Kglobs(b-2:b,b-2:b)+Kele(4:6,4:6);
    Kglobs(a-2:a,b-2:b)=Kglobs(a-2:a,b-2:b)+Kele(1:3,4:6);
    Kglobs(b-2:b,a-2:a)=Kglobs(b-2:b,a-2:a)+Kele(4:6,1:3);

end
end

