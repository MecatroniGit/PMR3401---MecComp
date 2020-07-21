clear all
clc

%Definição das variáveis necessárias na primeira tarefa
h=3;
d=5*h;
L=2*h;
H=8*h;
dx=h/6;
dy=dx;
V=100/3.6;
Comprimento=60;
ro=1.25;
gama=1.4;
w1=1.85; %fator de sobrerelaxação
erroe1=0.01;

%Inicializa a matriz de psi
psi=zeros(([(H/dy)+1 1+(2*d+L)/dx]));

%Define as coordenadas necessárias para os gráficos
[X,Y]=meshgrid(dx*((1:length(psi(1,:)))-1),((length(psi(:,1)):-1:1)-1)*dx);
x=15:0.01:21;

%Exercício 1
%Item a
[psi,bordairregular,AB]=ItemA(psi,H,h,d,L,dx,dy,V,w1,erroe1);

%Plota os gráficos do Item a
figure(1)
surf(X,Y,psi)
title("Função corrente do escoamento ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")
zlabel("Função corrente do escoamento")
colorbar

figure(2)
contour(X,Y,psi,40)
title("Função corrente do escoamento ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")
colorbar

%Item b
[u,v]=ItemB(psi,H,h,d,L,dx,dy,V);

%Plota os gráficos do Item b
figure(3)
quiver(X,Y,u,-v)
title("Velocidade do fluxo ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")

figure(4)
surf(X,Y,sqrt(u.^2+v.^2))
title("Módulo da velocidade ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")
zlabel("Módulo da velocidade (m/s)")
colorbar

%Item c
pressaovetor=ItemC(ro,gama,u,v);
 
%Plota o gráfico do Item c
figure(5)
surf(X,Y,pressaovetor)
title("Variação de pressão ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")
zlabel("Variação de pressão")
colorbar

%Item d
pressaoarrumada=pressao_telhado(psi,H,h,d,L,dx,dy,pressaovetor,bordairregular,AB);

%Plota os gráficos do Item d
figure(6)
plot(pressaoarrumada(2:end,1),pressaoarrumada(2:end,3))
title("Variação de pressão ao longo da borda")
xlabel("Eixo x (m)")
ylabel("Variação de pressão")

% figure(7)
% plot(x,sqrt(9-(x-18).^2)+3,'--'),grid
% hold on
% scatter(pressaoarrumada(2:end,1),pressaoarrumada(2:end,2),40,pressaoarrumada(2:end,3))
% hold on
% scatter(pontosoriginais(2:end,1),pontosoriginais(2:end,2),40,pressaoarrumada(2:end,3),'filled')
% title("Variação de pressão ao longo da borda")
% xlabel('Eixo x (m)')
% ylabel('Eixo y (m)')
% colorbar

%Item e
%[FR,Areas,Normais,Forca]=ItemE(pressaoarrumada,Comprimento);

%Plota os gráficos do Item e
% figure(8)
% plot(x,sqrt(9-(x-18).^2)+3,'--'),grid
% hold on
% scatter(pressaoarrumada(2:end,1),pressaoarrumada(2:end,2),40,pressaoarrumada(2:end,3),'filled')
% xlabel('Eixo x (m)')
% ylabel('Eixo y (m)')
% colorbar
% hold on
% quiver(pressaoarrumada(2:end,1),pressaoarrumada(2:end,2),Normais(:,1),Normais(:,2))
% title("Vetores normais ao longo da borda")
% 
% figure(9)
% plot(x,sqrt(9-(x-18).^2)+3,'--'),grid
% hold on
% scatter(pressaoarrumada(2:end,1),pressaoarrumada(2:end,2),40,pressaoarrumada(2:end,3),'filled')
% title("Força no telhado")
% xlabel('Eixo x (m)')
% ylabel('Eixo y (m)')
% colorbar
% hold on
% quiver(pressaoarrumada(2:end,1),pressaoarrumada(2:end,2),-Forca.*Normais(:,1),-Forca.*Normais(:,2))
% 
% display('A força resultante sobre o telhado, em Newtons, é')
% FR=FR

%Exercício 2
%Define as variáveis necessárias para a segunda tarefa que ainda não foram
%definidas
k=0.026;
cp=1002;
Tdentro=40;
Tfora=20;
w2=1; %Ajustado do enunciado pois para 1.15 não convergia
T=zeros(([(H/dy)+1 1+(2*d+L)/dx]));
T(:,:)=20;

%Item a
T=distribuicao_temp(T,k,ro,cp,Tdentro,Tfora,dx,dy,H,h,L,d,V,w2,u,v,erroe1);

%Plota os gráficos do Item a da segunda tarefa
figure(10)
surf(X,Y,T)
title("Temperatura ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")
zlabel("Temperatura (ºC)")
colorbar

figure(11)
contourf(X,Y,T+273.15)
title("Temperatura ao longo do domínio")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")
colorbar

%Item b
[calor,gradientesarrumado,Normais]=Item2b(T,bordairregular,AB,psi,dx,dy,H,h,L,d,Areas,Normais,k);

%Plota o gráfico do Item b da segunda tarefa
figure(12)
a=zeros([1 301]);
a(:,:)=15;
c=zeros([1 301]);
c(:,:)=21;
y=sqrt(9-(x-18).^2)+3;
z=[a x c];
w=[0:0.01:3 y 0:0.01:3];
plot(z,w,'--'),grid
hold on
quiver(gradientesarrumado(2:end,1),gradientesarrumado(2:end,2),gradientesarrumado(2:end,3).*Normais(:,1),gradientesarrumado(2:end,4).*Normais(:,2))
title("Vetores da transferência de calor ao longo do silo")
xlabel("Eixo x (m)")
ylabel("Eixo y (m)")

display('O calor total retirado do galpão, em W, é ')
calor=calor
