%Note: for files, be sure to remove any other text (ex. the name)
%Similar to gc_content... but different!

file1 = "/Users/lucywu/drylab/codon_harmony/onexoutput.fasta";
%replace the filepaths with chosen fasta file

%code is a little redundant below... but that can be easily fixed later

%matrixDisp(file1);
freq(file1);

function freq(input)
    file = fopen(input,'r');
    sequence = fscanf(file,'%c');

    fclose(file);
    
    sequence = regexprep(sequence,'\s','');
    sequence
    
    frequency = zeros(2, 20); %shouldn't really hard code this but...
    i = 2;
    f = 1; %Keep track of what part of matrix to append to

    while i < length(sequence)
        %Variable for count
        c = 1; 
        
        %If the previous is the same, increment
        while sequence(i) == 'C' || sequence(i) == 'G'
            c = c+1;
            i = i+1;
            c
        end
        %Otherwise, add to the matrix
        %Append the count
        f
        matrix(1, f) = c;

        %Append the letter
        matrix(2, f) = sequence(i);

        %Increment
        f = f+1;
        i = i+1;
    end
    frequency
    %Output the maximum consecutive Gs and Cs
    max(frequency)

end

function matrixDisp(input)
    file = fopen(input,'r');
    sequence = fscanf(file,'%c');

    fclose(file);
    
    sequence = regexprep(sequence,'\s','');
    
    matrix = zeros(1,length(sequence));
    for i =1:length(sequence)
        if sequence(i) == 'C' || sequence(i) == 'G'
            matrix(1,i) = 255;
        else
            matrix(1,i) = 0;
        end
    end
%{
%This was definitely not necessary...
    for i=2:(length(sequence)-1)
        if sequence(i-1) == 'C' || sequence(i-1) == 'G'
            if sequence(i+1) == 'C' || sequence(i+1) == 'G'
                %Surrounded on both sides by C and/or G
                matrix(1,i) = 2;
            else
                %Only preceding is a C or G
                matrix(1,i) = 1;
            end
        elseif sequence(i+1) == 'C' || sequence(i+1) == 'G'
            %Only successive is C or G
            matrix(1,i) = 1;
        else
            %Surrounding is not a C or G
            matrix(1,i) = 0;
        end
    end
%}
    matrix
    densityMap = reshape(matrix,1,[]);       
    figure;
    imshow(densityMap);
    customColorMap=jet(256);
    imagesc(densityMap);
    colormap(customColorMap);
    colorbar;

end