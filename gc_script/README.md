Motivation:
We are creating this script primarily to pass the IDT complexity test for ordering gene fragments: https://www.idtdna.com/site/order/gblockentry. Properly designing our primers is crucial to the success of our project.

One of the biggest challenges with ordering a silk protein sequence is that repeats and GC content become a big factor in an increasing complexity score. This is because a high GC nucleotide content can cause the templates to fold into more complex secondary structures during assembly and inhibit expression levels. Also, greater GC content contributes to higher melting temperatures.

1. To address length, we have used the codon optimization tools, such as the codon_harmony script, to eliminate unnecessary repeats and optimize the sequence for length.
2. However, GC content is still inflating our complexity score. Determining where to cut our sequence can help, since the proportion of G and C nucleotides in that segment can be minimized. Wetlab would then assemble the pieces together to form a whole sequence. 

Solution: From our literature review, there is no tool so far that accounts for this, possibly because this is not usually an issue. Thus, we will create a tool that optimizes the GC content of the segments.

Parameters:
1. Ideally, we should limit our segments to 16 pieces, since that is the max our wetlab team can handle assembling together. Increasing this parameter will make it more difficult to have a reasonable assembly.
2. GC content is calculated by taking the number of G and C nucleotides in the segment, and dividing it by the total number of nucleotides in the segment. Ideally, it should be kept between 65% and 40% in a segment.
3. A single segment may NOT have more than 3000 bp or less than 125 bp since that is the max and min set by IDT
4. We are assuming that your sequence would benefit from cutting: there is a tipping point where it would not be possible to have all pieces under 65% anyways, so we are assuming the sequence has a reasonably high and fixable GC content, at ~75%.

Goals:
1. Cut the sequence in such a way that all pieces are composed of less than 65% GC content. We should also keep GC content above 40%, otherwise GC content will be too low.
2. Minimize the number of pieces in the sequence to aid wetlab assembly. If possible, keep it under 16 pieces.
3. Arriving at a solution. Does not necessarily mean we are trying EVERY option, but with the power we have, we should get something.
Bonus: Minimizing the cost: if there is another solution that has a lower cost, list that as a better solution than another. Generally, having more short fragments is more expensive. https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments

Main Challenges:
1. Sequences are long. Like thousands of nucleotides long. This means there are an uncountable number of ways to cut up the sequence. While computing is fast, we have a limited amount of memory. So, we need to think about the most optimal way to arrive at a solution the quickest.
2. Storing sequences and trying something else?
3. Mathematically, there MUST be a point where even cutting can't save you. At what point would that be? Our program should tell you automatically based on the initial GC percentage to avoid wasting time and resources.

Strategies:
1. Piece by Piece: Start at the beginnning and add Gs and Cs to the piece until we hit %64. Repeat until there are no nucleotides left.
- Cons: good chance to still exceed 16 pieces. Also, this is definitely not the most efficient way.
2. Binary Search: Cut the sequence in half and measure GC content. If either fails, then cut that piece in half.
- Pros: much more efficient
- Cons: not the most fine method, may end up with pieces that are not optimal
3. Nile River: Locate the longest length of GC (or homopolyerimic run of G or C) and cut the sequence in half there. Rinse and repeat

Features:
1. Importing a fasta file with output (possibly not in fasta format)
2. Generating a random sequence to test
3. Incorporate other IDT parameters (GC near ends, etc.)

Questions:
1. Does there need to be some kind of overlap? Or can wetlab handle none
2. How do we feel about using... C++ instead of Python for some malloc abilities? - Let's prototype in Python first and if it gets computationally expensive... then lets move to C++
3. At the end, are we able to cut the n-terminus?