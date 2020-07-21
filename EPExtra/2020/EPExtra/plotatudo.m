function plotatudo(Vetor,parametro1,parametro2,extra1,extra2,...
    listanos,classeelementos)
    % Fun��o respons�vel por plotar os valores obtidos
    global h
    
    Plot = zeros(52/h, 26/h);
    for i=1:length(listanos)
      Plot(1+(listanos(i,2)+26)/h,1+listanos(i,1)/h)=Vetor(i);
    end
    
    Plotajeitado = [fliplr(Plot) Plot];
    Plotajeitado(:,length(Plotajeitado)/2)=[];     
    [X,Y]=meshgrid(-26:h:26, -26:h:26);
    
    figure();
    surf(X,Y,Plotajeitado);        
    xlabel('Eixo x (m)');
    ylabel('Eixo y (m)');  
    colorbar
    if length(extra1)==1
        zlabel('V (kV)');  
        title('Potencial el�trico V');
        salvar=strcat('Vh',num2str(h),'.jpg');
        saveas(gcf,salvar)
    else
        zlabel('Az (Wb/m)');  
        title('Potencial vetor magn�tico Az');
        salvar=strcat('Az',num2str(h),'.jpg');
        saveas(gcf,salvar)
    end

    
    figure();
    contourf(X, Y,Plotajeitado); 
    [xmedio,ymedio]=baricentro(classeelementos);
    hold on
    quiver([-fliplr(xmedio) xmedio],[fliplr(ymedio) ymedio],...
    [-fliplr(parametro1) parametro1],[fliplr(parametro2) parametro2],...
    0.1,'k');
    colorbar 
    xlabel('Eixo x (m)');
    ylabel('Eixo y (m)');
    if length(extra1)==1
        title('Potencial el�trico V e campo el�trico E');
        salvar=strcat('VEh',num2str(h),'.jpg');
        saveas(gcf,salvar)
    else
        title('Potencial vetor magn�tico Az e campo magn�tico B')
        salvar=strcat('ABh',num2str(h),'.jpg');
        saveas(gcf,salvar)
    end
    axis image;
    
    if length(extra1)~=1
        figure();
        quiver([-fliplr(xmedio) xmedio],[fliplr(ymedio) ymedio],...
        [-fliplr(extra1) extra1],[fliplr(extra2) extra2],0.1,'k');
        title('Campo magn�tico H');
        xlabel('Eixo x (m)');
        ylabel('EIxo y (m)');
        axis image;
        salvar=strcat('Hh',num2str(h),'.jpg');
        saveas(gcf,salvar)
    end
    
end

