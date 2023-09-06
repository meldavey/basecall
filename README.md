Base caller example
Mel Davey, 2023

example usage:

performs phase correction, calls bases, shows plots:
python3 ./basecall.py -ie 0.09 -cf 0.065 -dr 0.02 -plots

performs a grid-search over a hard-coded range of cf/ie/dr params, then performs base calling on best cf/ie/dr params
python3 ./basecall.py -grid -plots

model.py - models the physical system
I was playing around with dark bases, but I think they would incorporate at a different rate, and not be UV cleave dependent so it needs work
not included here yet

basecaller.py - the basic app to perform base calling as follows:
1. loads a dataset of measured signals for each dye at the start of each cycle
2. instantiates a model to calculate expected signals based on phase params
3. and applies the model to predict signals using a what-if strategy at each cycle, trying each base and seeing which was best
4. best base is set in the template (initially blank), and a multi-pass approach can make better use of the model and partial
   template to account for carry-forward effects

jmrdata1.json - a json file containing the dye intensity measurements at each cycle
This was a fun dataset, and requires a little background.  Jonathan Rothberg had tweeted a puzzle of sorts, which caught my attension, given that I had previously solved this sort of thing at Ion Torrent.  The DNA being sequenced is a small group of identical copies per spot/well/etc.  This is not a single molecule sequencing.  So, over time, when the polymerase imcorporates a base or the DNA strands expose the next nucleotide on the strand, not all polymerases advance, some stay behind, some can incorporate multiple bases, etc.  So over time, the strands get "out of phase".  It's context-dependent and is complex to solve.  So I've developed a simple prediction model that can predict the signal measured based on the current phasing of polymerases, and one of 4 possible next bases.  The system simply makes predictions, and picks the best one.  The other part of the trick is that the model cannot initially predict the signal generated from polymerases advancing beyond the current base prediction.  So a multi-pass approach is used, where the first pass gets a few bases correct where it can account for signal loss and polymerases that fail to advance, but that gives us enough information to predict the future on the second and subsequent passes.  This turns out to be pretty fast (well, when written correctly in C).  So, back to the dataset...   I grabbed it from his tweet, converted it into the json file seen here, and applied the phase and signal loss corrections, and got a really good prediction of the bases modeled in the original case.  My solution was much more accurate than any other attempt.   I've since given the team my code (I'm retaining copywrite on it for now), which started with what I posted here.

future work:
1. fix the dark nuc model
2. expand into a dynamic programming style, so that a single incorrect base call won't impact the entire system downstream of that call,
   basically track the top-N best calls, or keep all 4 base calls for say the last 4 cycles, dropping the lowest scoring paths as
   each new cycle is processed.  Might be better to just track the top-10 DNA template possibilities and prune at each new cycle

