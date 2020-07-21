function Derivadas = geraDerivadas(L1,L2,R,m1,m2,miIz,Vel,L2eixo,F1,F2,g)
%Definindo os termos necessários para as derivadas
%i será um vetor-coluna, com os valores necessários de Theta1, Theta2,
%Alpha1 e Alpha2, respectivamente.
%t é o instante de tempo que devemos considerar

    A0=@(t,i) L1^2*L2*R*(m2*cos(2*i(1)-2*i(2))-2*m1-m2);
    A1=@(t,i) L1^2*L2*R*m2*sin(2*i(1)-2*i(2));
    A2=@(t,i) 2*L1*L2^2*R*m2*sin(i(1)-i(2));
    A3=@(t,i) -2*L2*miIz*Vel;
    A4=@(t,i) -2*L1*miIz*Vel*cos(i(1)-i(2));
    A5=@(t,i) -R*L1*(L2eixo*F2*sin(i(1)-2*i(2))+2*sin(i(1))*(F1*L2+L2eixo*F2*0.5));

    B0=@(t,i) L2^2*R*m2;
    B1=@(t,i) -L1*L2*R*m2*cos(i(1)-i(2));
    B2=@(t,i) L1*L2*R*m2*sin(i(1)-i(2));
    B3=@(t,i) -miIz*Vel;
    B4=@(t,i) L2eixo*sin(i(2))*R*F2;

    %Derivadas - será usado o vetor de variáveis no lugar de i na hora de
    %integrar
    Theta1_Ponto=@(t,i) i(3);
    Theta2_Ponto=@(t,i) i(4);
    Alpha1_Ponto=@(t,i) (A1(t,i)*Theta1_Ponto(t,i)^2+A2(t,i)*Theta2_Ponto(t,i)^2+A3(t,i)*i(3)+A4(t,i)*i(4)+A5(t,i))/A0(t,i);
    Alpha2_Ponto=@(t,i) (B1(t,i)*Alpha1_Ponto(t,i)+B2(t,i)*Theta1_Ponto(t,i)^2+B3(t,i)*Theta2_Ponto(t,i)+B4(t,i))/B0(t,i);

    Derivadas= @(t,i) [Theta1_Ponto(t,i);Theta2_Ponto(t,i);Alpha1_Ponto(t,i);Alpha2_Ponto(t,i)];
end

