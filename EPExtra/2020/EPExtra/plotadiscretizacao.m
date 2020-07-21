function plotadiscretizacao(listanos,classeelementos)
    % Função responsável por plotar os elementos e os respectivos nós para
    % confirmar que está correto
    global h
    figure();
    x=-26:0.01:26;
    plot(x,sqrt(26^2-x.^2),'k','LineWidth',2)
    hold on
    plot(x,-sqrt(26^2-x.^2),'k','LineWidth',2)
    for i=1:length(classeelementos)
        hold on
        plot([classeelementos(i).no1posx classeelementos(i).no2posx],...
            [classeelementos(i).no1posy classeelementos(i).no2posy],'b')    
        plot(-[classeelementos(i).no1posx classeelementos(i).no2posx],...
            [classeelementos(i).no1posy classeelementos(i).no2posy],'b')
        plot([classeelementos(i).no2posx classeelementos(i).no3posx],...
            [classeelementos(i).no2posy classeelementos(i).no3posy],'b')
        plot(-[classeelementos(i).no2posx classeelementos(i).no3posx],...
            [classeelementos(i).no2posy classeelementos(i).no3posy],'b')
        plot([classeelementos(i).no3posx classeelementos(i).no1posx],...
            [classeelementos(i).no3posy classeelementos(i).no1posy],'b')
        plot(-[classeelementos(i).no3posx classeelementos(i).no1posx],...
            [classeelementos(i).no3posy classeelementos(i).no1posy],'b')
    end
    scatter(listanos(:,1),listanos(:,2),7.5*h,'r','filled')
    hold on
    scatter(-listanos(:,1),listanos(:,2),7.5*h,'r','filled')
    xlabel('Eixo x (m)')
    ylabel('Eixo y (m)')
    titulo=strcat('Discretização da malha h=',num2str(h));
    title(titulo)
    axis image
    salvar=strcat('discretizacaoh',num2str(h),'.jpg');
    saveas(gcf,salvar)
    
end

