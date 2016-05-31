# Laser Albatross

*Author:* Martín Basterrechea

## Done

* Support for multiple input and output formats
* It works for both nucleotide and amino acid sequences.
* Print rejected blocks
* Print accepted and rejected block ranges
* Support for gaps
* Support for HTML output
* HTML output shows which blocks have been removed
* "All, partial and none" gap support:
  * *None:* Gaps are considered non-conserved positions and are taken out.
  * *Partial:* Gaps are considered just another allowed character.
  * *All:* Gaps positions are not removed, and they don't count towards the conserved or non conserved calculation of the position. This means that the required amount of identical characters for a position is lower.
* Info-score info in HTML output




## Gap support

An example for the difference of partial and all using the standard settings ("conserved" is (0.5 * # of sequences) + 1  and "highly conserved" is 0.85 * # of seqs.) in this case 5 and 7, respectively.
In this alignment:


    ---
    ---
    VWQ
    AWL
    VWQ
    VWQ
    AWQ
    ---


Position 1 has 3 gaps. With *partial* setting they are considered another character. No characters are above the *conserved* threshold so this position is *non conserved*. With *all* settings, the threshold is lower ((0.5 * 5 ) + 1), 3.5, but no other character is above this threshold so this is *not conserved*.

Position 2 has 3 gaps. With *partial* setting the most common character is **W**, and it is just in the threshold (5), so this position is considered *conserved*. With *all* settings, the gaps are not considered a character and the threshold goes down (as they don't count towards the threshold calculation), so the new thresholds for this position are 3.5 and 4.25. **W** is present in 5 sequences so this position is considered *Highly conserved*.

Position 3 has 3 gaps. It would be considered *non conserved* in *partial* setting and *conserved* in *all* setting.




## Todo
* Refactor code
* Maybe have both alignments in the html scroll at the same time.
* Maybe integrate with bash pipes? check emboss
* Add support for ambiguous characters
* Add polydot chart?
* Clean params/args stuff
* Automatic score threshold calculation
* Score looks really interesting, maybe it can be added a chemical group score and then average.



## Useful info

[Information content](http://www.lecb.ncifcrf.gov/~toms/paper/primer/)
