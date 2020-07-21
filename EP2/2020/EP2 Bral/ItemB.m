function [u,v]=ItemB(psi,H,h,d,L,dx,dy,V)
%Resolve o Item B da Parte 1 do EP2

%Inicializa as matrizes das derivadas
dpsidx=zeros(([(H/dy)+1 1+(2*d+L)/dx]));
dpsidy=zeros(([(H/dy)+1 1+(2*d+L)/dx]));

for j=1:length(psi(:,1))
    for i=1:(length(psi(1,:))+1)/2
        if i==1 || j==1
            %Confere se est� em alguma das laterais
            if j~=1        
                dpsidx(j,i)=0;
                if j~=length(psi(:,1))
                    dpsidy(j,i)=(psi(j+1,i)-psi(j-1,i))/(2*dx);
                else
                    %Est� na quina superior esquerda
                    dpsidy(j,i)=(psi(j,i)-psi(j-1,i))/dx; 
                end
            elseif i~=1
                dpsidx(j,i)=(psi(j,i+1)-psi(j,i-1))/(2*dx);
                dpsidy(j,i)=-V;
            else
                dpsidx(j,i)=0;
                dpsidy(j,i)=-V;                
            end
        elseif j==length(psi(:,1))
            dpsidx(j,i)=(psi(j,i+1)-psi(j,i-1))/(2*dx);
            dpsidy(j,i)=(psi(j,i)-psi(j-1,i))/dx;
        elseif i>=((1/dx)*d+1) && j>=(1/dx)*(H-h)+1
            %Confere se est� dentro do silo (parte retangular)
            dpsidx(j,i)=0;
            dpsidy(j,i)=0;
        elseif j<(1/dx)*(H-h)+1 && j>=(1/dx)*(H-L/2-h)+1 && i>=(1/dx)*d+1
            %Estamos no ret�ngulo no qual o semicirculo do silo est�
            %inscrito. O cume, em L/2, ter� uma dist�ncia dx com rela��o ao
            %�ltimo ponto acima dele sempre, pois ser�o escolhidos valores
            %de dx e dy para tal, assim como a dist�ncia para a lateral do
            %silo
            if (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j-1))^2<=(0.5*L)^2
                %Estamos dentro do semicirculo, e portanto, devemos definir
                %esses pontos como zero
                %Checar direitinho se ta pegando tudo que era pra ta aqui
                %dentro
                if (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j-1))^2==(0.5*L)^2
                    dpsidx(j,i)=0;
                    dpsidy(j,i)=(psi(j-2,i)-4*psi(j-1,i)+3*psi(j,i))/...
                        (2*dx);   
                else
                    dpsidx(j,i)=0;
                    dpsidy(j,i)=0;
                end
            elseif (dx*(i)-d-0.5*L)^2+(H-h-dx*(j-1))^2<(0.5*L)^2 || ...
                    (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j))^2<(0.5*L)^2
                %Checa se o pr�ximo ponto vai cair no (tanto em x quanto em
                %y) cai dentro do c�rculo. Se for cair, temos que usar o MDF 
                %para bordas irregulares para calcular a itera��o neles.
                a=(d+L/2-dx*(i-1)-sqrt((L/2)^2-(H-h-dx*(j-1)).^2))/dx;
                b=(H-h-dx*(j-1)-sqrt((L/2)^2-(d+0.5*L-dx*(i-1))^2))/dx;
                if a<1 && b<1
                    dpsidx(j,i)=(psi(j,i+1)-a^2*psi(j,i-1)-(1-a^2)*...
                        psi(j,i))/((a+a^2)*dx);
                    dpsidy(j,i)=(psi(j+1,i)-b^2*psi(j-1,i)-(1-b^2)*...
                        psi(j,i))/((b+b^2)*dx);
                elseif a<1
                    dpsidx(j,i)=(psi(j,i+1)-a^2*psi(j,i-1)-(1-a^2)*...
                        psi(j,i))/((a+a^2)*dx);
                    dpsidy(j,i)=(psi(j+1,i)-psi(j,i))/(dx);
                elseif b<1
                    dpsidx(j,i)=(psi(j,i+1)-psi(j,i))/dx;
                    dpsidy(j,i)=(psi(j+1,i)-b^2*psi(j-1,i)-(1-b^2)*...
                        psi(j,i))/((b+b^2)*dx);
                end
            else
                dpsidx(j,i)=(psi(j,i+1)-psi(j,i-1))/(2*dx);
                dpsidy(j,i)=(psi(j+1,i)-psi(j-1,i))/(2*dx);      
            end
        else
            dpsidx(j,i)=(psi(j,i+1)-psi(j,i-1))/(2*dx);
            dpsidy(j,i)=(psi(j+1,i)-psi(j-1,i))/(2*dx);            
        end
        if i<=(length(psi(1,:)))/2 
            dpsidx(j,length(psi(1,:))+1-i)=-dpsidx(j,i);
            %Inverte o sinal pois em vez de o vento estar desviando do silo
            %est� voltando a ocupar o espa�o no qual havia o silo.
            dpsidy(j,length(psi(1,:))+1-i)=dpsidy(j,i);
        elseif i>=length(psi(1,:))/2
            dpsidx(j,length(psi(1,:))+1-i)=dpsidx(j,i);
            dpsidy(j,length(psi(1,:))+1-i)=dpsidy(j,i);  
        end  
    end
end

%Salva os vetores das derivadas como os componentes das velocidades,
%fazendo os ajustes de sinais necess�rios pela orienta��o do sistema
u=-dpsidy;
v=dpsidx;
end

