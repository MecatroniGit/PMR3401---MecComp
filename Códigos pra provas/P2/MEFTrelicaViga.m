%Esse c�digo funciona pra no m�ximo 2 elementos. Em nenhuma prova o Em�lio
%cobrou mais do que 2, ent�o duvido que d� a louca nesse.

%Insere aqui os do problema (se for s� 1 elemento, s� precisa modificar os
%elementos 1
E1=10;
I1=0.1;
A1=2;
L1=2;
ro1=1;

E2=10;
I2=0.1;
A2=2;
L2=2;
ro2=1;

%Insere os �ngulos aqui, em radianos
theta1=0;
theta2=90*pi/180;

%Insere aqui o n�mero de elementos
Num_elementos=2;

%Insere 1 se for treli�a ou 2 se for viga
Elemento=2;

%Insere abaixo a matriz de carregamentos globais

%Insere aqui a an�lise que voce quer fazer
%1- Harm�nica 2- Modal 3- Est�tica (pra treli�a)
analise=2;
%Se for harm�nica, descomenta e coloca o valor da frequencia
%w=2;

%Insere aqui os carregamentos (se for modal, comenta isso)


if Elemento==1
    %Cria as matrizes de cada elemento
    Kele=@(E,A,L) E*A/L *[1 0 -1 0; 0 0 0 0; -1 0 1 0; 0 0 0 0];
    Kele=@(ro,A,L) ro*A*L/6 *[2 0 1 0; 0 0 0 0; 1 0 2 0; 0 0 0 0];
    %Cria as matrizes T com os thetas
    T1=[cos(theta1) sin(theta1) 0 0; -sin(theta2) cos(theta2) 0 0;...
        0 0 cos(theta1) sin (theta2); 0 0 -sin(theta2) cos(theta2)];
    T2=[cos(theta2) sin(theta2) 0 0; -sin(theta2) cos(theta2) 0 0;...
        0 0 cos(theta2) sin (theta2); 0 0 -sin(theta2) cos(theta2)];
else
    %Cria as matrizes de cada elemento
    Kele=@(E,I,L) E*I*L^3 * [12 6*L -12 6*L; 6*L 4*L^2 -6*L 2*L^2; -12 -6*L 12 -6*L; 6*L 2*L^2 -6*L 4*L^2];
    Mele=@(ro,A,L) ro*A*L/420 *[156 22*L 54 -13*L; 22*L 4*L^2 13*L -3*L^2;54 13*L 156 -22*L;-13*L -3*L^2 -22*L 4*L^2];
end

%Monta matrizes de cada um dos elementos e a global
if Num_elementos==1
    K=T1'*Kele(E1,I1,L1)*T1; %Inserir aqui o valor de E,I,L do enunciado
    M=T1'*Mele(ro1,A1,L1)*T1; %Inserir ro,A,L do enunciado
else
    K1=K(E1,I1,L1);
    K2=K(E2,I2,L2);
    Kglobal=zeros(6);
    Kglobal(1:4,1:4)=K1;
    Kglobal(3:6,3:6)=Kglobal(3:6,3:6)+K2;

    M1=M(ro1,A1,L1);
    M2=M(ro2,A2,L2);
    Mglobal=zeros(6);
    Mglobal(1:4,1:4)=M1;
    Mglobal(3:6,3:6)=Mglobal(3:6,3:6)+M2;
end




