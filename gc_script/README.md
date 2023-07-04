Motivation:
We are creating this script primarily to pass the IDT complexity test for ordering gene fragments: https://www.idtdna.com/site/order/gblockentry. Properly designing our primers is crucial to the success of our project.

One of the biggest challenges with ordering a silk protein sequence is that repeats and GC content become a big factor in an increasing complexity score. This is because a high GC nucleotide content can cause the templates to fold into more complex secondary structures during assembly and inhibit expression levels. Also, greater GC content contributes to higher melting temperatures.

1. To address length, we have used the codon optimization tools, such as the codon_harmony script, to eliminate unnecessary repeats and optimize the sequence for length.
2. However, GC content is still inflating our complexity score. Determining where to cut our sequence can help, since the proportion of G and C nucleotides in that segment can be minimized. Wetlab would then assemble the pieces together to form a whole sequence. 

Solution: From our literature review, there is no tool so far that accounts for this, possibly because this is not usually an issue. Thus, we will create a tool that optimizes the GC content of the segments.

Parameters:
1. Ideally, we should limit our segments to 16 pieces, since that is the max our wetlab team can handle assembling together. Increasing this parameter will make it more difficult to have a reasonable assembly.
2. GC content is calculated by taking the number of G and C nucleotides in the segment, and dividing it by the total number of nucleotides in the segment.

Goals:
1. Cut the sequence in such a way that all pieces are composed of less than 65% GC content. We should also keep GC content above 40%, otherwise GC content will be too low.
2. Minimize the number of pieces in the sequence to aid wetlab assembly. If possible, keep it under 16 pieces.
3. Arriving at a solution. Does not necessarily mean we are trying EVERY option, but with the power we have, we should get something.
Bonus: Minimizing the cost: if there is another solution that has a lower cost, list that as a better solution than another. Generally, having more short fragments is more expensive. https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments

Main Challenges:
1. Sequences are long. Like thousands of nucleotides long. This means there are an uncountable number of ways to cut up the sequence. While computing is fast, we have a limited amount of memory and ways to try. So, we need to think about the most optimal way to arrive at a solution the quickest.
2. Storing sequences and trying something else?

Stratgies:
1. 

Questions:
1. Does there need to be some kind of overlap? Or can wetlab handle none
2. How do we feel about using... C instead of Python for some malloc abilities?