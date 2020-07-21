function Forca= geraF(nos_F1,nos_F2,t,tamanho,modo)
    %Fun��o para gerar o vetor das for�as pra an�lise transiente e
    %harm�nica
    global t1 t2 F D
    
    Forca=zeros(tamanho,1);
    if modo==1
        Forca(3*nos_F1-1)=2*F*sin(2*pi*t);
        if t>=t1 && t<=t2
            Forca(3*nos_F2-2)=5*D;
        end
    else
        Forca(3*nos_F1-1)=2*F;
    end
end

