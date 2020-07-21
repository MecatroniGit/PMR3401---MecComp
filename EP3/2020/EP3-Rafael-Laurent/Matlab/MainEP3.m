close all

% Definição das variáveis
global F D t1 t2 L1 L2 L3 v d1i d1e d2i d2e E ro Elementos 
global alfa beta gama beta2
F=8000;
D=2000;
t1=2;
t2=8;
L1=2;
L2=3;
L3=4;
v=0.29;
d1i=0.072;
d1e=0.08;
d2i=0.09;
d2e=0.1;
E=210e9;
ro=7650;
Elementos=4;

alfa=0.3;
beta=0.03;
gama=1/2;
beta2=1/4;


tic
%Pré-processamento - Definição de nós e barras
nos_originais=[0 0;4 0;0.5 3;3.5 3;1 6;3 6;1 8;3 8;1 10;3 10;1 12;3 12;1 14;3 14];
barras=criabarras(nos_originais);

%Cria matrizes por elemento
[K,M,Malternativa]=criamatrizes;

%Cria os elementos e atualiza os nós
[nos_elementos,elementos]=criaelementos(barras);

%Cria as matrizes globais
[Kglobs,Mglobs,Mglobs2,Lele]=matrizesglobais(elementos,nos_elementos,K,M,Malternativa);

toc

%Análise modal
tic
[A,w1,CC]=modal(nos_elementos,elementos,Kglobs,Mglobs,Lele);
toc
%Análise modal com a matriz 'lumped'
tic
[B,w2,CC]=modal(nos_elementos,elementos,Kglobs,Mglobs2,Lele);
toc

%Análise transiente
tic
[MCC,KCC,nos_F1,nos_F2] = transiente(nos_elementos,CC,Mglobs,Kglobs);
toc
%Análise transiente com a matriz 'lumped'
tic
[MCC2,KCC,nos_F1,nos_F2] = transiente(nos_elementos,CC,Mglobs2,Kglobs);
toc

%Análise harmônica
tic
harmonica(nos_elementos,MCC,KCC,nos_F1,nos_F2)
toc
%Análise harmônica com a matriz 'lumped'
tic
harmonica(nos_elementos,MCC2,KCC,nos_F1,nos_F2)
toc






