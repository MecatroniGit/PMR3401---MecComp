%Passo 1 - Discretizar o domínio
global nos elementos dTdn Tinf hbarra Lbarra ki dTdntipo nosCC TCC

%Inserir aqui as coordenadas dos nos no espaço
nos=[0 0; 1 0; 0 1; 1 1];

%Nós que conheco a temperatura
%Temperatura conhecida
nosCC=[0 0 0 1];
TCC=[0 0 0 100];

%Inserir aqui, na ordem antihoraria, os nos dos elementos
%Caso seja triangulo e nao quadrilatero, por um 0 no fim
%Caso seja uma barra, por 2 zeros no fim
%Se for quadrado, começa do extremo inferior esquerdo
elementos=[1 2 3 0; 2 4 0 0];

%Inserir aqui os k de cada elemento
ki=[2 3];

%Insere se tem carregamento distribuido no elemento ou não
Carregamento=[0];
FCar=2;

%Insere indices das condições de contorno por contorno
%Se for triangular ou barra, nas posicoes extras poe dTdntipo=-1
%dT/dn=cte -> dTdntipo=1
%kdT/dn=-h(T-Tinf) -> dTdntipo=2
%dT/dn=-a*T -> dTdntipo=2 também
%dT/dn=desconhecido -> dTdntipo=3
%dT/dn=interface -> dTdntipo=4
%Primeiro elemento é contorno do nó 1 pro 2 e assim por diante
dTdntipo=[1 4 3 1; 1 3 4 -1];
%Aqui coloca o valor de dTdn, caso seja constante.
%Coloca o valor de a, caso seja do tipo 3. Coloca o valor de h caso seja
%tipo 2
%Caso seja tipo -1, 3 ou 4, poe 1000
dTdn=[2 1000 1000 0; 0 1000 1000 -1];

%Se não tiver, poe Inf pra saber onde da matriz por a convecção
Tinf=0;

%Convecção em barra - caso tenha - se não for barra só põe 0 e foda se
hbarra=[0 0];
Lbarra=1;

Kglobal=zeros(length(nos(:,1)));
Kglobal2=zeros(length(nos(:,1)));
Fglobal=zeros(length(nos(:,1)),length(find(elementos~=0)));
bc=zeros(3,2);
nosconveccao=[];
Kglobsadd=[];

for i=1:length(elementos(:,1))
    elemento=find(elementos(i,:)==0);
    if (isempty(elemento)) %quadrilatero
        L=nos(elementos(i,2),1)-nos(elementos(i,1),1)
        W=nos(elementos(i,4),2)-nos(elementos(i,1),2)
        Kele1=ki(i)*L/(6*W)*[2 1 -1 -2;1 2 -2 -1; -1 -2 2 1; -2 -1 1 2];
        Kele2=ki(i)*W/(6*L)*[2 -2 -1 1; -2 2 1 -1; -1 1 2 -2; 1 -1 -2 2];
        disp('A primeira matriz (L/W) é')
        Kele1
        disp('A segunda matriz (W/L) é')
        Kele2
        disp('E a matriz final é')
        Kele=Kele1+Kele2
        [Fele,nosconveccao,Kglobsadd]=F(Kele,i,Kglobal,nosconveccao,Kglobsadd)
    elseif length(elemento)==1 %triangulo
        bc(1,1)=nos(elementos(i,2),2)-nos(elementos(i,3),2);
        bc(2,1)=nos(elementos(i,3),2)-nos(elementos(i,1),2);
        bc(3,1)=nos(elementos(i,1),2)-nos(elementos(i,2),2);
        bc(1,2)=nos(elementos(i,3),1)-nos(elementos(i,2),1);
        bc(2,2)=nos(elementos(i,1),1)-nos(elementos(i,3),1);
        bc(3,2)=nos(elementos(i,2),1)-nos(elementos(i,1),1);
        disp('Os valores de b e c são:')
        bc
        disp('E a area calculada é:')
        Ae=0.5*(bc(1,1)*bc(2,2)-bc(2,1)*bc(1,2))
        Kele=zeros(3,3);
        Kele(1,:)=[dot(bc(1,:),bc(1,:)) dot(bc(1,:),bc(2,:)) dot(bc(1,:),bc(3,:))];
        Kele(2:3,1)=Kele(1,2:3)';
        Kele(2,2:3)=[dot(bc(2,:),bc(2,:)) dot(bc(2,:),bc(3,:))];
        Kele(3,2)=Kele(2,3);
        Kele(3,3)=dot(bc(3,:),bc(3,:));
        disp('A matriz K do triangulo é')
        Kele
        disp('Depois de multiplicar por k/4Ae fica:')
        Kele=ki(i)/(4*Ae)*Kele
        [Fele,nosconveccao,Kglobsadd]=F(Kele,i,Kglobal,nosconveccao,Kglobsadd);
    else %barra pq o emilio eh um lixo fudido
        disp('A matriz K da barra é:')
        Seila=[1 -1;-1 1]
        disp('Multiplicando por k/L fica')
        Kele1=ki(i)/Lbarra*[1 -1; -1 1]
        disp('A matriz H da barra é:')
        seila=[1 0.5;0.5 1];
        disp('E a matriz H multiplicando por hL/3 é:')
        Kele2=hbarra(i)*(Lbarra/3)*[1 0.5; 0.5 1]
        disp('E a K final é')
        Kele=Kele1-Kele2
        Fele=0;
    end
    disp('A F desse elemento é:')
    Fele
    Kglobal=incrementaKglobal(Kglobal,Kele,i);
    if length(elemento)<=1
        Fglobal=incrementaFglobal(Fglobal,Fele,i);
    end
end

disp('A Kglobal fica:')
Kpreglobs=Kglobal
disp('E a Fglobal fica:')
Fglobal
disp('E a Fglobal ajeitada fica:')
Fglobs=AjeitaFglobal(Fglobal)

if length(nosconveccao)~=0
    Kglobal2=mudaKglobal(Kglobal,nosconveccao,Kglobsadd);
    disp('Apos a conveccao, a matriz global fica:')
    Kposglobs=Kglobal2
else
    Kglobal2=Kglobal;
end
[Fglobs,Fglobs2,Kglobal,Kglobal2]=Montasistema(Kglobal,Kglobal2,Fglobs);

disp('A matriz Kglobal, tirando os GL nulos, fica:')
Kglobal2

disp('A matriz F, tirando os GL nulo, fica:')
Fglobs2

Resolvesistema(Fglobs2,Kglobal2);

function Resolvesistema(Fglobs2,Kglobal2)
    global nosCC
    Kmini=zeros(2,2);
    Fmini=zeros(2,1);
    manter=find(nosCC==0);
    Kmini(1,1)=Kglobal2(manter(1),manter(1));
    Kmini(1,2)=Kglobal2(manter(1),manter(2));
    Kmini(2,1)=Kglobal2(manter(2),manter(1));
    Kmini(2,2)=Kglobal2(manter(2),manter(2));
    Fmini(1)=Fglobs2(manter(1));
    Fmini(2)=Fglobs2(manter(2));
    disp('A matriz K, diminuida pra 2x2, fica:')
    Kmini
    disp('A matriz F, diminuida pra 2x2, fica:')
    Fmini
    disp('E a solução final é:')
    A=Kmini\Fmini
end

function [Fglobs,Fglobs2,Kglobal,Kglobal2]=Montasistema(Kglobal,Kglobal2,Fglobs)
    global nosCC TCC
    Fglobs2=Fglobs
    for i=1:length(Kglobal)
        if nosCC(i)==1
            Fglobs=Fglobs-Kglobal(:,i)*TCC(i);
            Kglobal(i,:)=0;
            Kglobal(:,i)=0;
            Kglobal(i,i)=1;
            Fglobs2=Fglobs2-Kglobal2(:,i)*TCC(i);
            Kglobal2(i,:)=0;
            Kglobal2(:,i)=0;
            Kglobal2(i,i)=1;            
        end
    end
    Fglobs(find(nosCC==1))=TCC(find(nosCC==1));
    Fglobs2(find(nosCC==1))=TCC(find(nosCC==1));
end

function Kglobal2=mudaKglobal(Kglobal,nosconveccao,Kglobsadd)
    Kglobal2=Kglobal;
    Kglobal2(nosconveccao(1),nosconveccao(1))=...
    Kglobal2(nosconveccao(1),nosconveccao(1))+Kglobsadd(1,1);
    Kglobal2(nosconveccao(1),nosconveccao(2))=...
    Kglobal2(nosconveccao(1),nosconveccao(2))+Kglobsadd(1,2);
    Kglobal2(nosconveccao(2),nosconveccao(1))=...
    Kglobal2(nosconveccao(2),nosconveccao(1))+Kglobsadd(2,1);
    Kglobal2(nosconveccao(2),nosconveccao(2))=...
    Kglobal2(nosconveccao(2),nosconveccao(2))+Kglobsadd(2,2);

end

function Fglobs=AjeitaFglobal(Fglobal)
    Fglobs=zeros(length(Fglobal(:,1)),1);
    for i=1:length(Fglobal(:,1))
        for j=1:length(Fglobal(1,:))
            if Fglobal(i,j)~=9999 && Fglobal(i,j)~=Inf && (~isnan(Fglobal(i,j)))
                Fglobs(i)=Fglobs(i)+Fglobal(i,j);
            end
        end
    end
end

function Fglobal=incrementaFglobal(Fglobal,Fele,i)
    global elementos
    for j=1:length(elementos)
        forma=length(find(elementos(i,:)==0));
        forma1=length(find(elementos(1,:)~=0));
        for k=1:length(elementos)-forma
            if elementos(i,j)~=0 && elementos(i,k)~=0
                m=elementos(i,j);
                n=k;
                if i==1
                    Fglobal(m,n)=Fglobal(m,n)+Fele(j,k);
                else
                    Fglobal(m,n+forma1)=Fglobal(m,n+forma1)+Fele(j,k);
                end
            end
        end
    end
end

function Kglobal=incrementaKglobal(Kglobal,Kele,i)
    global elementos
    for j=1:length(elementos)
        for k=1:length(elementos)
            if elementos(i,j)~=0 && elementos(i,k)~=0
                m=elementos(i,j);
                n=elementos(i,k);
                Kglobal(m,n)=Kglobal(m,n)+Kele(j,k);
            end
        end
    end
end

%Deixar -1000 na matriz F de onde tem que cortar
%Substituir quando for por na global
function [Fele,nosconveccao,Kglobsadd]=F(Kele,i,Kglobal,nosconveccao,Kglobsadd)
    global dTdn Tinf hbarra Lbarra elementos nos ki dTdntipo
    Fele=ones(length(Kele));
    for j=1:length(Kele) %Ns
        for k=1:length(Kele) %Contornos
            zerei=0;
            forma=find(elementos(i,:)==0);
            if k+1<=length(elementos(i,:))-length(forma)
                if elementos(i,j)~=elementos(i,k) && elementos(i,j)~=elementos(i,k+1)
                    %Tem que zerar porque é Nj no contorno Lik
                    Fele(j,k)=0;
                    zerei=1;
                end
            elseif elementos(i,j)~=elementos(i,k) && elementos(i,j)~=elementos(i,1)
                %Tem que zerar porque é Nj no contorno Lik
                Fele(j,k)=0;
                zerei=1;
            end
            if zerei==0
                if dTdntipo(i,k)==1 
                    %dT/dn=cte -> dTdntipo=1
                    cte=dTdn(i,k);
                    if k+1<=length(elementos(i,:))-length(forma)
                        if nos(elementos(i,k),1)==nos(elementos(i,k+1),1)
                            dij=abs(nos(elementos(i,k),2)-nos(elementos(i,k+1),2));
                        else
                            dij=abs(nos(elementos(i,k),1)-nos(elementos(i,k+1),1));
                        end
                    else
                        if nos(elementos(i,k),1)==nos(elementos(i,1),1)
                            dij=abs(nos(elementos(i,k),2)-nos(elementos(i,1),2));
                        else
                            dij=abs(nos(elementos(i,k),1)-nos(elementos(i,1),1));
                        end                       
                    end
                    calc=ki(i)*cte*(dij)/2; 
                    Fele(j,k)=calc;
                elseif dTdntipo(i,k)==2
                    %dT/dn=-h(T-Tinf) -> dTdntipo=2
                    %dT/dn=-a*T -> dTdntipo=2
                    if k+1<=length(elementos(i,:))-length(forma)
                        if length(nosconveccao)==0
                            nosconveccao=[elementos(i,k) elementos(i,k+1)];
                        end
                        if nos(elementos(i,k),1)==nos(elementos(i,k+1),1)
                            Fele(j,k)=Tinf*dTdn(i,k)*0.5*abs(nos(elementos(i,k),2)-nos(elementos(i,k+1),2));
                            if length(Kglobsadd)==0
                                Kglobsadd=dTdn(i,k)*abs(nos(elementos(i,k),2)-nos(elementos(i,k+1),2))*[1/3 1/6; 1/6 1/3];
                            end
                        else
                            Fele(j,k)=Tinf*dTdn(i,k)*0.5*abs(nos(elementos(i,k),1)-nos(elementos(i,k+1),1));  
                            if length(Kglobsadd)==0
                                Kglobsadd=dTdn(i,k)*abs(nos(elementos(i,k),1)-nos(elementos(i,k+1),1))*[1/3 1/6; 1/6 1/3];
                            end                        
                        end
                    else
                        if length(nosconveccao)==0
                            nosconveccao=[elementos(i,k) elementos(i,1)];
                        end
                        if nos(elementos(i,k),1)==nos(elementos(i,1),1)
                            Fele(j,k)=Tinf*dTdn(i,k)*0.5*abs(nos(elementos(i,k),2)-nos(elementos(i,1),2));
                            if length(Kglobsadd)==0
                                Kglobsadd=dTdn(i,k)*abs(nos(elementos(i,k),2)-nos(elementos(i,1),2))*[1/3 1/6; 1/6 1/3];
                            end                        
                        else
                            Fele(j,k)=Tinf*dTdn(i,k)*0.5*abs(nos(elementos(i,k),1)-nos(elementos(i,1),1));
                            if length(Kglobsadd)==0
                                Kglobsadd=dTdn(i,k)*abs(nos(elementos(i,k),1)-nos(elementos(i,1),1))*[1/3 1/6; 1/6 1/3];
                            end                        
                        end   
                    end
                elseif dTdntipo(i,k)==3
                    %dT/dn=desconhecido -> dTdntipo=3
                    Fele(j,k)=NaN;
                elseif dTdntipo(i,k)==4
                    %dT/dn=interface -> dTdntipo=4 
                    Fele(j,k)=9999;
                end
            end
        end
    end

end