function [psi,bordairregular,AB]=ItemA(psi,H,h,d,L,dx,dy,V,w1,erroe1)
%Resolve o Item A da Parte 1 do EP2

bordairregular=[0 0];
contador=0;
AB=[0 0];

contiteracoes=0;
itera=1;

while itera==1
    errovelho=0; %Zera o erro para comparar a cada iteração 
    erronovo=0;
    for j=1:(length(psi(:,1))-1)
        %Não precisa ir até o final porque a última linha possui condição de
        %Dirichlet, fixa em psi=0, então não itera nela
        for i=1:(length(psi(1,:))+1)/2
            %Vai somente até o meio porque é simétrico o problema
            psivelho=psi(j,i);
            if i==1
                %Devemos aplicar as condições de contorno onde dpsi/dx=0
                if j==1
                    %Condição de contorno para a quina superior esquerda, que
                    %possui duas condições de Neumann
                    psi(j,i)=(psi(j,i+1)+psi(j+1,i)+V*dx)/2;
                else
                    psi(j,i)=(2*psi(j,i+1)+psi(j-1,i)+psi(j+1,i))/4;
                end
            elseif j==1
                %Devemos aplicar as condições de contorno onde dpsi/dy=-V
                psi(j,i)=(2*psi(j+1,i)+psi(j,i-1)+psi(j,i+1)+2*V*dx)/4;
            elseif i>=((1/dx)*d+1) && j>=(1/dx)*(H-h)+1
                %Confere se está dentro do silo (parte retangular)
                psi(j,i)=0;
            elseif j<(1/dx)*(H-h)+1 && j>=(1/dx)*(H-L/2-h)+1 && i>=(1/dx)*d+1
                %Estamos no retângulo no qual o semicirculo do silo está
                %inscrito. O cume, em L/2, terá uma distância dx com relação ao
                %último ponto acima dele sempre, pois serão escolhidos valores
                %de dx e dy para tal, assim como a distância para a lateral do
                %silo
                if (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j-1))^2<=(0.5*L)^2
                    %Estamos dentro do semicirculo, e portanto, devemos definir
                    %esses pontos como zero
                    %Checar direitinho se ta pegando tudo que era pra ta aqui
                    %dentro
                    psi(j,i)=0;
                    if contiteracoes==0 && (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j-1))^2==(0.5*L)^2
                        bordairregular(end+1,:)=[j,i];
                        %Salva os índices dos valores que serão necessários
                        %e as distâncias A e B deles para futuro uso nos
                        %itens 1D, 1E e 2B
                        AB(end+1,:)=[2,2]; %distancia do ponto a esquerda
                        contador=contador+1;
                    end
                elseif (dx*(i)-d-0.5*L)^2+(H-h-dx*(j-1))^2<(0.5*L)^2 || ...
                        (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j))^2<(0.5*L)^2
                    %Checa se o próximo ponto vai cair no (tanto em x quanto em
                    %y) cai dentro do círculo. Se for cair, temos que usar o MDF 
                    %para bordas irregulares para calcular a iteração neles.
                    a=(d+L/2-dx*(i-1)-sqrt((L/2)^2-(H-h-dx*(j-1)).^2))/dx;
                    b=(H-h-dx*(j-1)-sqrt((L/2)^2-(d+0.5*L-dx*(i-1))^2))/dx;
                    if contiteracoes==0
                        %Salva os índices dos valores que serão necessários
                        %e as distâncias A e B deles para futuro uso nos
                        %itens 1D, 1E e 2B
                        bordairregular(end+1,:)=[j,i];                        
                        bordairregular(end+1,:)=[j,length(psi(1,:))+1-i];
                        contador=contador+1;
                        AB(end+1,:)=[a,b]; %distancia do ponto a esquerda
                        AB(end+1,:)=[a,b]; %distancia do ponto a direita
                    end                    
                    if a<1 && b<1
                        %Ponto está mais perto que dx e que dy do circulo
                        psi(j,i)=(psi(j,i+1)/(a*(a+1))+psi(j,i-1)/(1+a)+psi(j+1,i)/(b*(b+1))...
                            +psi(j-1,i)/(1+b))/(1/a+1/b);
                    elseif a<1
                        %Ponto está a uma distância maior ou igual que dy
                        %mas menor que dx do círculo
                        psi(j,i)=(2*psi(j,i+1)/(a+a^2)+2*psi(j,i-1)/(1+a)+psi(j+1,i)...
                            +psi(j-1,i))/(2+2/a);
                    elseif b<1
                        %Ponto está a uma distância maior ou igual que dx
                        %mas menor que dy do círculo
                        psi(j,i)=(2*psi(j+1,i)/(b+b^2)+2*psi(j-1,i)/(1+b)+psi(j,i+1)...
                            +psi(j,i-1))/(2+2/b);
                    end
                else
                    psi(j,i)=(psi(j+1,i)+psi(j-1,i)+psi(j,i+1)+psi(j,i-1))/4;
                end
            else
               psi(j,i)=(psi(j+1,i)+psi(j-1,i)+psi(j,i+1)+psi(j,i-1))/4;
            end
            psi(j,i)=w1*psi(j,i)+(1-w1)*psivelho;
            if i<=(length(psi(1,:)))/2-1
                psi(j,length(psi(1,:))+1-i)=psi(j,i);
            elseif i>=length(psi(1,:))/2
                psi(j,length(psi(1,:))+2-i)=psi(j,i-1);
            end
            erronovo=abs(psi(j,i)-psivelho);
            if erronovo>errovelho
                %Armazena o maior valor do erro (em módulo)
                errovelho=erronovo;
            end
        end
    end
    if errovelho<=erroe1
        %Caso o máximo erro calculado seja menor que o aceitavel, encerra o
        %loop
        itera=0;
    end
    contiteracoes=contiteracoes+1; 
    %Conta quantas iterações foram necessárias pra convergir
end

display('O total de iterações necessário para esse passo foi')
contiteracoes=contiteracoes
end

