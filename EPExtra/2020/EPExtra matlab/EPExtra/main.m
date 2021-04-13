clear all
close all

%% Definição dos valores usados
global mi1 mi2 sigma1 sigma2 h raio Imaximo RC fonte1 fonte2 Vmaximo borda

fonte1=[4,14];
fonte2=[6,10];
h=2;   
raio=26; 
Imaximo=200;
RC=0.02;
Vmaximo=500;
borda=1;
mi1=1.2566*10^-6;
mi2=2*1.2566*10^-6;
sigma1=10^-10;
sigma2=10^-2;

%% Pré-processamento - Criação dos elementos e dos nós
classeelementos=criaelementos;
[listanos,classeelementos]=crianos(classeelementos);

%% Pré-processamento - Estabelecimento das Condições de Contorno
[pontosborda,pontosfonte,pontoscarro]=geraCCs(listanos);

%% Loop que resolve tanto para V e obtém E 
% quanto resolve para Az e obtém B e H
for parametro=1:2
    %% Pré-processamento-Criação das matrizes de rigidez e de carregamento
    K=criaglobal(listanos,classeelementos,parametro);
    [F,K]=novaF(listanos,pontosborda,pontosfonte,pontoscarro,K,parametro);

    %% Resolução do sistemma 
    Solucao=linsolve(K, F);
    
    %% Pós-processamento
    [Vetorx,Vetory]=solveEB(classeelementos,Solucao,parametro);
    
    if parametro==2
        Hx=zeros(length(Vetorx),1);
        Hy=zeros(length(Vetory),1);
        for i = 1:length(classeelementos)
            if classeelementos(i).no3posy>=0 &&...
            classeelementos(i).no2posy>=0 && classeelementos(i).no1posy>=0
                Hx(i) = Vetorx(i)/mi1;
                Hy(i) = Vetory(i)/mi1;
            else
                Hx(i)= Vetorx(i)/mi2;
                Hy(i)= Vetory(i)/mi2;
            end
        end
    else
        Hx=1;
        Hy=1;
    end
    
    %% Plotagem do resultado
    plotatudo(Solucao,Vetorx,Vetory,Hx,Hy,listanos,classeelementos);
end








