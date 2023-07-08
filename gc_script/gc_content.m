file1 = "C:\Users\ninja\Downloads\output_onex_05relax.fasta";
file2 = "C:\Users\ninja\Downloads\output_onex.fasta";
%replace the filepaths with chosen fasta file

matrixDisp(file1);
matrixDisp(file2);




function matrixDisp(input)
    file = fopen(input,'r');
    sequence = fscanf(file,'%c');

    fclose(file);
    
    sequence = regexprep(sequence,'\s','');
    
    matrix = zeros(1,length(sequence));
   

    for i =1:length(sequence)
        if sequence(i) == 'C' || sequence(i) == 'G'
            row = 1;
            col = i;
            matrix(row,col) = 255;
        end
    end
    densityMap = reshape(matrix,1,[]);       
    figure;imshow(densityMap);
    customColorMap=jet(256);
    imagesc(densityMap);
    colormap(customColorMap);
    colorbar;
    

end