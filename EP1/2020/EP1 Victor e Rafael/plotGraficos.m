function plotGraficos(titulo,tempo,Variaveis,Arquiva_Derivadas)
% plotGraficos = Plota o gráfico desejado, na forma de subplot para
% facilitar a visualização, com legenda e título.

    figure('Name',titulo,'NumberTitle','off')
    title(titulo);
    %TAPE MUDA PRA SGTITLE ESSA PORRA

    subplot(3,2,1);
    plot(tempo,Variaveis(1,:));
    grid on;
    xlabel('Tempo (s)')
    ylabel('$\theta_{1} (rad)$','interpreter','latex','Rotation',0)

    subplot(3,2,2);
    plot(tempo,Variaveis(2,:));
    grid on;
    xlabel('tempo (s)')
    ylabel('$\theta_{2} (rad)$','interpreter','latex','Rotation',0)

    subplot(3,2,3);
    plot(tempo,Variaveis(3,:));
    grid on;
    xlabel('tempo (s)')
    ylabel('$\dot{\theta_{1}} (rad)$','interpreter','latex','Rotation',0)

    subplot(3,2,4);
    plot(tempo,Variaveis(4,:));
    grid on;
    xlabel('tempo (s)')
    ylabel('$\dot{\theta_{2}} (rad)$','interpreter','latex','Rotation',0)

    subplot(3,2,5);
    plot(tempo,Arquiva_Derivadas(3,:));
    grid on;
    xlabel('tempo (s)')
    ylabel('$\ddot{\theta_{1}} (rad)$','interpreter','latex','Rotation',0)

    subplot(3,2,6);
    plot(tempo,Arquiva_Derivadas(4,:));
    grid on;
    xlabel('tempo (s)')
    ylabel('$\ddot{\theta_{2}} (rad)$','interpreter','latex','Rotation',0)

end

