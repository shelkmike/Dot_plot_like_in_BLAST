## Introduction
The online BLAST (https://blast.ncbi.nlm.nih.gov) makes dot plots which look like this:  
&nbsp;  
<img src="http://mikeshelk.site/Diff/Files_for_GitHub/Dot_plot_like_in_BLAST/dot_plot_from_online_BLAST.png?" width="500">
&nbsp;  
&nbsp;  
Unfortunately, the *standalone* BLAST lacks the capability of making dot plots. I have made a tool named **Dot_plot_like_in_BLAST** specifically for this purpose. Basically, it is a standalone tool that aligns by BLAST two nucleotide sequences and makes a dot plot. The main differences from the online BLAST are:  
1) While the online BLAST always makes rectangular images with a 2:1 aspect ratio, the aspect ratio in Dot_plot_like_in_BLAST is equal to the ratio of sequences' lengths. In other words, Dot_plot_like_in_BLAST doesn't distort lengths.  
2) The online BLAST draws direct and reverse matches in one color. Dot_plot_like_in_BLAST draws direct matches in blue, and reverse matches in orange.  
3) Dot_plot_like_in_BLAST has several options for image customization. For example, it allows to change the font size and line width.  
4) The online BLAST make images in PNG. Dot_plot_like_in_BLAST makes images in PNG, SVG, and also it makes interactive files in HTML. In the interactive file you can magnify a region of interest and save an SVG file specifically for this region.  
&nbsp;  
This is how an image made by Dot_plot_like_in_BLAST looks. It was made for the same sequences as the image above.
&nbsp;  
&nbsp;  
<img src="http://mikeshelk.site/Diff/Files_for_GitHub/Dot_plot_like_in_BLAST/Dot_plot_like_in_BLAST.png?" width="500">
&nbsp;  
&nbsp;  
An example of interactive image manipulation (a region is magnified and saved as a separate SVG file):
<img src="http://mikeshelk.site/Diff/Files_for_GitHub/Dot_plot_like_in_BLAST/Dot_plot_like_in_BLAST__animation.gif?" width="500">

## How to install
Dot_plot_like_in_BLAST has two prerequisites:  
1) The path to BLAST+ executable files should be in $PATH . BLAST+ can be downloaded here https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/  
2) The Python library Plotly should be installed. You can install it, for example, using the following command:  
`pip3 install --upgrade --user plotly`  
&nbsp;  
Then, download the latest version of Dot_plot_like_in_BLAST from [Releases](https://github.com/shelkmike/Dot_plot_like_in_BLAST/releases) and extract the archive. That's all.  
&nbsp;  
Dot_plot_like_in_BLAST works on all operating systems.  

## How to use
Dot_plot_like_in_BLAST has only two mandatory options:
1. --file_with_the_first_sequence . This is the path to a FASTA file with the first sequence.  
2. --file_with_the_second_sequence . This is the path to a FASTA file with the second sequence.

<!-- Empty comment to end a list. Otherwise the list with three options will have numbers "3.", "4.", "5.". See https://stackoverflow.com/a/58506405 -->

A typical command of Dot_plot_like_in_BLAST looks like this:  
`python3 dot_plot_like_in_BLAST.py --file_with_the_first_sequence sequence_1.fasta --file_with_the_second_sequence sequence_2.fasta`  
&nbsp;  
Dot_plot_like_in_BLAST creates a folder with results. The three main files in this folder are: image.png, image.svg, image.html. They are, respectively, a raster image with the dot plot, a vector image with the dot plot, and an interactive image with the dot plot.  
&nbsp;  
Dot_plot_like_in_BLAST has many options that allow you to customize alignment and visualization. For example:  
1. --minimum_percent_identity . Sets the minimum percent identity for a match.  
2. --minimum_match_length . Sets the minimum match length.  
3. --font_size . Sets the font size in the image.  
&nbsp;  
And many others. To see the full list of options, run `python3 dot_plot_like_in_BLAST.py --help`

## Questions and answers
1) Where can I suggest features and report problems?  
To suggest features and report problems, please post in Issues (https://github.com/shelkmike/Dot_plot_like_in_BLAST/issues). 
2) Can Dot_plot_like_in_BLAST be used for very long sequences, like entire eukaryotic chromosomes?  
Yes, but the calculation can take a lot of time if you use the default parameters. You may need to increase --minimum_percent_identity and --minimum_match_length; decrease --number_of_best_matches_to_draw; set --blast_tool_to_use to "megablast".  
For example, with the following settings the calculation will be faster than with default settings: "--minimum_percent_identity 80 --minimum_match_length 1000 --blast_tool_to_use megablast".
3) How does Dot_plot_like_in_BLAST depict indels?  
Just as in the online BLAST, indels (gaps in the alignment) are depicted by gaps in the diagram.  
For example, here is an alignment (a part of a longer alignment) and a corresponding diagram:  
<img src="http://mikeshelk.site/Diff/Files_for_GitHub/Dot_plot_like_in_BLAST/alignment_with_a_gap_and_its_diagram.png?" width="500">  
