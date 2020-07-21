function T=Item2a(T,k,ro,cp,Tdentro,Tfora,dx,dy,H,h,L,d,V,w2,u,v,erroe1)
%Resolve o Item A da Parte 2 do EP2

itera=1;
contiteracoes=0;

while itera
    errovelho=0;
    erronovo=0;
    for j=1:(length(T(:,1)))
        for i=1:(length(T(1,:)))
            Tvelho=T(j,i);
            if i==1 || j==1 || i==length(T(1,:)) || j==length(T(:,1))
                %Condições de contorno do domínio
                if i==1
                    T(j,i)=Tfora;
                elseif j==1
                    if i~=length(T(1,:))
                        %Parte superior do domínio
                        %Só consideramos o caso de u>0 pois é o que
                        %acontece em todo o domínio
                        T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+2*T(j+1,i))...
                               +ro*cp*(u(j,i)*T(j,i-1)))/...
                                (4*k/dx+ro*cp*(u(j,i)));      
                    else
                        %Quina do lado superior direito
                        T(j,i)=0.5*(T(j,i-1)+T(j+1,i));
                    end
                elseif i==length(T(1,:))
                    %Parte direita do domínio
                    if j~=length(T(:,1))
                        %v=0, então vai cancelar independente do sinal
                        T(j,i)=(k/dx*(2*T(j,i-1)+T(j-1,i)+T(j+1,i))...
                            )/(4*k/dx);   
                    else
                        %Quina do lado inferior direito
                        T(j,i)=0.5*(T(j,i-1)+T(j-1,i));
                    end                    
                elseif j==length(T(:,1))
                    %Chão do domínio
                    if i<(1/dx)*d+1 || i>(1/dx)*(d+L)+1
                        %Fora do hangar
                        %Só consideramos o caso de u>0 pois é o que
                        %acontece em todo o domínio
                        T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+2*T(j-1,i))...
                               +ro*cp*(u(j,i)*T(j,i-1)))/...
                                (4*k/dx+ro*cp*(u(j,i)));       
                    else
                        T(j,i)=Tdentro;
                    end
                end
            elseif i>=((1/dx)*d+1) && j>=(1/dx)*(H-h)+1 && i<=(1/dx)*...
                    (d+L)+1
                %Confere se está dentro do silo (parte retangular)
                T(j,i)=Tdentro;
            elseif j<(1/dx)*(H-h)+1 && j>=(1/dx)*(H-L/2-h)+1 && i>=...
                    (1/dx)*d+1 && i<=(1/dx)*(d+L)+1
                %Estamos no retângulo no qual o semicirculo do silo está
                %inscrito. O cume, em L/2, terá uma distância dx com relação ao
                %último ponto acima dele sempre, pois serão escolhidos valores
                %de dx e dy para tal, assim como a distância para a lateral do
                %silo
                if (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j-1))^2<=(0.5*L)^2
                    %Estamos dentro do semicirculo, e portanto, devemos definir
                    %esses pontos como zero
                    T(j,i)=Tdentro;
                elseif (dx*(i)-d-0.5*L)^2+(H-h-dx*(j-1))^2<(0.5*L)^2 || ...
                        (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j))^2<(0.5*L)^2 || ...
                        (dx*(i-2)-d-0.5*L)^2+(H-h-dx*(j-1))^2<(0.5*L)^2                        
                    %Checa se o próximo ponto - ou anterior - vai cair 
                    %dentro do círculo(tanto em x quanto em y). Se for 
                    %cair, temos que usar o MDF para bordas irregulares 
                    %para calcular a iteração neles.
                    if (dx*(i)-d-0.5*L)^2+(H-h-dx*(j-1))^2<(0.5*L)^2  
                        %Proximo ponto em x é o hangar
                        if (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j))^2<(0.5*L)^2
                            %Proximo ponto em y é o hangar
                            a=(d+L/2-dx*(i-1)-sqrt((L/2)^2-...
                                (H-h-dx*(j-1)).^2))/dx;
                            b=(H-h-dx*(j-1)-sqrt((L/2)^2-...
                                (d+0.5*L-dx*(i-1))^2))/dx;
                            %Usando a aproximação aconselhada para a
                            %primeira derivada de T, pelos resultados
                            %anteriores sabemos que nesse caso, u>0 e v<0
                            %em todos os pontos. Assim, temos que
                            T(j,i)=(2*k/dx*(T(j,i+1)/(a+a^2)+T(j,i-1)/...
                                (1+a)+T(j+1,i)/(b^2+b)+T(j-1,i)/(b+1))-...
                                ro*cp*(-u(j,i)*T(j,i-1)+v(j,i)*...
                                T(j+1,i)/b))/(2*k/dx*(1/a+1/b)+ro*cp*...
                                (u(j,i)-v(j,i)/b));
                        else
                            a=(d+L/2-dx*(i-1)-sqrt((L/2)^2-...
                                (H-h-dx*(j-1)).^2))/dx;
                            T(j,i)=(k/dx*(2*T(j,i+1)/(a+a^2)+2*T(j,i-1)/...
                                (a+1)+T(j+1,i)+T(j-1,i))-ro*cp*(-u(j,i)...
                                *T(j,i-1)+v(j,i)*T(j+1,i)))/(k/dx*(2+2/a)...
                                +ro*cp*(u(j,i)-v(j,i)));
                        end
                    elseif (dx*(i-2)-d-0.5*L)^2+(H-h-dx*(j-1))^2<(0.5*L)^2 
                        %Ponto anterior em x é o hangar
                        if (dx*(i-1)-d-0.5*L)^2+(H-h-dx*(j))^2<(0.5*L)^2
                            %Proximo ponto em y é o hangar
                            a=(dx*(i-1)-(d+L/2)-sqrt((L/2)^2-...
                                (H-h-dx*(j-1)).^2))/dx;                            
                            b=(H-h-dx*(j-1)-sqrt((L/2)^2-...
                                (d+0.5*L-dx*(i-1))^2))/dx;
                            T(j,i)=(2*k/dx*(T(j,i+1)/(a+1)+T(j,i-1)/...
                            (a+a^2)+T(j+1,i)/(b+1)+T(j-1,i)/(b^2+b))-...
                            ro*cp*(-u(j,i)*T(j,i-1)/a-v(j,i)*T(j-1,i)))/...
                            (2*k/dx*(1/a+1/b)+ro*cp*(u(j,i)/a+v(j,i)));
                        else
                            a=(dx*(i-1)-d-L/2-sqrt((L/2)^2-...
                                (H-h-dx*(j-1)).^2))/dx;
                            T(j,i)=(k/dx*(2*T(j,i+1)/(a+1)+2*T(j,i-1)/...
                            (a+a^2)+T(j+1,i)+T(j-1,i))-ro*cp*(-u(j,i)*...
                            T(j,i-1)/a-v(j,i)*T(j-1,i)))/(k/dx*(2+2/a)+...
                            ro*cp*(u(j,i)/a+v(j,i)));                            
                        end
                    else
                        %Proximo ponto em y (somente) é o hangar
                        if v(j,i)<=0
                            b=(H-h-dx*(j-1)-sqrt((L/2)^2-...
                                (d+0.5*L-dx*(i-1))^2))/dx;
                            T(j,i)=(k/dx*(2*T(j+1,i)/(b+b^2)+...
                            2*T(j-1,i)/(b+1)+T(j,i+1)+T(j,i-1))-...
                            ro*cp*(-u(j,i)*T(j,i-1)+v(j,i)*T(j+1,i)/b))/...
                            (k/dx*(2+2/b)+ro*cp*(u(j,i)-v(j,i)/b));
                        else
                            T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+T(j-1,i)+...
                            T(j+1,i))+ro*cp*(u(j,i)*T(j,i-1)+v(j,i)*...
                            T(j-1,i)))/(4*k/dx+ro*cp*(u(j,i)+v(j,i)));                            
                        end
                    end
                else
                    %Pontos internos sem problemas no domínio
                    %Consideramos u>0 sempre pois é o que acontece em todo
                    %o domínio
                    if v(j,i)>=0
                        T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+T(j-1,i)+T(j+1,i))...
                           +ro*cp*(u(j,i)*T(j,i-1)+v(j,i)*T(j-1,i)))/...
                            (4*k/dx+ro*cp*(u(j,i)+v(j,i)));
                    else
                        T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+T(j-1,i)+T(j+1,i))...
                           +ro*cp*(u(j,i)*T(j,i-1)-v(j,i)*T(j+1,i)))/...
                            (4*k/dx+ro*cp*(u(j,i)-v(j,i)));                    
                    end
                end      
            else
                %Pontos internos sem problemas no domínio
                %u>0 em todo o domínio
                if v(j,i)>=0
                    T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+T(j-1,i)+T(j+1,i))...
                       +ro*cp*(u(j,i)*T(j,i-1)+v(j,i)*T(j-1,i)))/...
                        (4*k/dx+ro*cp*(u(j,i)+v(j,i)));
                else
                    T(j,i)=(k/dx*(T(j,i+1)+T(j,i-1)+T(j-1,i)+T(j+1,i))...
                       +ro*cp*(u(j,i)*T(j,i-1)-v(j,i)*T(j+1,i)))/...
                        (4*k/dx+ro*cp*(u(j,i)-v(j,i))); 
                end
            end
            T(j,i)=w2*T(j,i)+(1-w2)*Tvelho;
            erronovo=abs(T(j,i)-Tvelho);
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
contiteracoes
end

