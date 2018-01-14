# Deduper

### Part 1
Write up a strategy for writing a reference-based PCR duplicate removal tool. Be sure to include an example input and output SAM file.

[Here](https://htmlpreview.github.io/) is a link to view the `deduper_part1.html`.

### Part 2

Wrote comments on Adrian's, Zach's, and Brandon's code.

### Part 3

Write the Python script for Deduper.

What's in this repository:


File | Comment
------ | -----------------
`version1_do_not_grade` | Directory containing old version of deduper script, Jupyter notebook with testing notes, and test files (for storage, not for grading).
`Deduper_lecture.pdf` | Leslie's lecture from 9 October 2017.
`STL96.txt` | UMI reference list for `Dataset1.sam`, `Dataset2.sam`, `Dataset3.sam`, and `paired_end.sam`.
`claridge_deduper.py` | Python script to run deduper.
`deduper.srun` | Slurm script to sort the SAM files and run deduper script, but my script runs quickly enough on the command line.
`deduper_part1.Rmd` | Deduper Part 1 pseudocode.
`deduper_part1.html` | Deduper Part 1 pseudocode.