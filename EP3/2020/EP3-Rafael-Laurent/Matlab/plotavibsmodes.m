function plotavibsmodes(A,nos_elementos,elementos,Lele)
    %Função para realizar a plotagem dos modos de vibração
    global d2e
    
    A=[zeros(1,6); zeros(1,6); A];
    posicao=find(nos_elementos(:,1)==4 & nos_elementos(:,2)==0);
    A=-[A(1:3*posicao-2,:); zeros(1,6); A(3*posicao-1:end,:)];

    xvet=zeros(length(elementos),21);
    yvet=zeros(length(elementos),21);

    nosnovos=nos_elementos;
    for i=1:length(A(1,:))
        ipx=20; %Numero de pontos dentro do elemento usados para interpolar
        figure()
        plotatorre(nos_elementos,elementos,3,i)
        for j=1:length(elementos) 
            %Percorre todos os elementos pra encontrar seus autovetores
            xf=[]; 
            %Vetor para armazenar os pontos dentro do elemento 
            %que serão aproximados por Hermite/Lagrange
            Df=[]; %Vetor para armazenar os 
            %autovetores associados só aquele elemento
            achax1=find(nos_elementos(:,1)==elementos(j,1));
            achay1=find(nos_elementos(achax1,2)==elementos(j,2));
            achax2=find(nos_elementos(:,1)==elementos(j,3));
            achay2=find(nos_elementos(achax2,2)==elementos(j,4));
            a=3*achax1(achay1);
            b=3*achax2(achay2);
            Autovet(1:3)=A(a-2:a,i); 
            %Salva as coordenadas associadas a esses pontos pra esse modo
            Autovet(4:6)=A(b-2:b,i);
            nosnovos(a/3,1)=nos_elementos(a/3,1)+Autovet(1);
            nosnovos(a/3,2)=nos_elementos(a/3,2)+Autovet(2);
            nosnovos(b/3,1)=nos_elementos(b/3,1)+Autovet(4);
            nosnovos(b/3,2)=nos_elementos(b/3,2)+Autovet(5);
            cosseno=(elementos(j,3)-elementos(j,1))/Lele(j);
            seno=(elementos(j,4)-elementos(j,2))/Lele(j);       
            for ip=1:ipx+1
                x=Lele(j)*(ip-1)/ipx; 
                %pega o ponto x do elemento discretizado pra interpolar
                L=Lele(j);
                %Lagrange
                N(1)=(L-x)/L;
                N(2)=x/L;
                %Hermite
                N(3)=1-3*x^2/L^2+2*x^3/L^3;
                N(4)=x-2*x^2/L+x^3/L^2;
                N(5)=3*x^2/L^2-2*x^3/L^3;
                N(6)=-x^2/L+x^3/L^2;
                Df(ip)=N(3)*Autovet(2)+N(4)*Autovet(3)+...
                    N(5)*Autovet(5)+N(6)*Autovet(6);
                xf(ip)=x+N(1)*Autovet(1)+N(2)*Autovet(4);
            end
            if i==1
            pos=find(nos_elementos(:,1)==elementos(j,1) &...
                nos_elementos(:,2)==elementos(j,2));
            xvet(j,:)=nosnovos(pos,1)+xf*cosseno-Df*seno;
            yvet(j,:)=nosnovos(pos,2)+Df*cosseno+xf*seno;
            else
            xvet(j,:)=elementos(j,1)+xf*cosseno-Df*seno;
            yvet(j,:)=elementos(j,2)+Df*cosseno+xf*seno;
            end
            if elementos(j,6)==d2e
                plot(xvet(j,:),yvet(j,:),'r','LineWidth',3)
                hold on
            else
                plot(xvet(j,:),yvet(j,:),'r','LineWidth',1)
                hold on
            end
        end
        hold on
        if i==1
            title("Tower's 1st Vibration Mode")
        elseif i==2
            title("Tower's 2nd Vibration Mode")
        elseif i==3
            title("Tower's 3rd Vibration Mode")
        elseif i>3
            titulo=strcat("Tower's ",num2str(i),"th Vibration Mode");
            title(titulo)
        end 
        xlabel('X axis');
        ylabel('Y axis');
        axis image
    end

end

