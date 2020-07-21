function plotatorre(nos, barras,titulo,modo)
    %Função para plotar a torre sem deformações
    global d2e
    if titulo~=3
       figure()

        for i = 1:length(barras)
            if barras(i,6) == d2e
                plot([barras(i,1) barras(i,3)],...
                    [barras(i,2) barras(i,4)],'LineWidth',3);
            else
                plot([barras(i,1) barras(i,3)],...
                    [barras(i,2) barras(i,4)],'LineWidth',1);
            end
            hold on
        end
        for i = 1:length(nos)
           scatter(nos(i,1),nos(i,2),'filled');
           hold on
        end

        xlabel('X axis');
        ylabel('Y axis');
        if titulo==1
            title('Tower original geometry');
        else
            title('Tower geometry with the division of the elements')
        end
        axis image
        hold off 
        if titulo==1
        saveas(gcf,'imgbarras.jpg')
        else
            saveas(gcf,'imgelementos.jpg')
        end
    else
        for i = 1:length(barras)
            if barras(i,6) == d2e
                plot([barras(i,1) barras(i,3)],...
                    [barras(i,2) barras(i,4)],'--k','LineWidth',2);
            else
                plot([barras(i,1) barras(i,3)],...
                    [barras(i,2) barras(i,4)],'--k','LineWidth',1);
            end
            hold on    
        end
    end
end