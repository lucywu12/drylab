%Note: for files, be sure to remove any other text (ex. the name)
%Similar to gc_content... but different!

file1 = "/Users/lucywu/drylab/codon_harmony/onexoutput.fasta";
%replace the filepaths with chosen fasta file

%freq(file1);
percent(file1);

function percent(input)

    %maybe rename some variables, there is something wrong with the
    %counting again

    file = fopen(input,'r');
    seq = fscanf(file,'%c');

    fclose(file);
    
    seq = regexprep(seq,'\s','');

    sequence = zeros(1, numel(seq)); %preallocate size of vector

    %create a sequence where GC is 1 and AT is 0
    for s=1:numel(seq)
        if seq(s) == 'G' || seq(s) == 'C'
            sequence(s) = 1;
        else
            sequence(s) = 0;
        end
    end

    ratios = zeros(1, length(sequence) - 125 + 1); %set up a place to store values
    
    %iterate through all 125 base pair frames
    
    count = 0;
    for i = 1:(length(sequence) - 125 + 1)
        count = 0;
        for frame = i:(i + 125 - 1)
            if sequence(frame) == 1 %if G or C
                
                count = count + 1; %increment count
            end
            ratios(i) = count; %add to our matrix
        end
    end

    ratios
    figure;
    bar(ratios)
end



function freq(input)

    file = fopen(input,'r');
    seq = fscanf(file,'%c');

    fclose(file);
    
    seq = regexprep(seq,'\s','');

    sequence = zeros(1, numel(seq)); %preallocate size of vector

    %create a sequence where GC is 1 and AT is 0
    for s=1:numel(seq)
        if seq(s) == 'G' || seq(s) == 'C'
            sequence(s) = 1;
        else
            sequence(s) = 0;
        end
    end

    
    frequency = zeros(3, 1); %start off
    i = 1;
    f = 1; %Keep track of what part of matrix to append to

    while i < length(sequence)
        length(sequence)
        %Variable for count
        c = 1; 
        
        %If the previous is the same, increment
        while (sequence(i) == sequence (i + 1)) && (i < length(sequence) - 1)
            c = c+1;
            i = i+1;

        end
        %Otherwise, add to the matrix
        %Append the count
        frequency(1, f) = c;

        %Append the letter (GC = 1, AT = 0)
        frequency(2, f) = sequence(i);

        %Append the base position
        frequency(3, f) = i - c + 1;

        %Increment
        f = f+1;
        i = i+1;
        
    end
    %For the last nucleotide
    if sequence(end) == sequence(end - 1)
        frequency(1, f-1) = c + 1;
    else 
        frequency(1, f) = 1
        frequency(2, f) = sequence(end)
        frequency(3, f) = i
    end

    %Output the maximum consecutive Gs and Cs
    %max(frequency)
    
    %frequency(1,:) %first row

    %plot the results
    figure; %generate a new window
    bar(frequency(1,:));
end

%Output: Count, GC/AT, Position
