This repo has everything you need to get some data on the variants within genes.

# Running The Program Variant_Finder.py (in terminal):
### General:
  * python2 Variant_Finder.py -h --help
  
  * usage: Variant_Finder.py [-h] [-s START] [-e END] in_File ex_File col_name {transcripts,genes}


### Positional Arguments:

  * in_File               File containing ID's you would like to find variants for.
                                               
  * ex_File               File exported from program with all the significant variants for provided IDs.
                                           
  * col_name              Name of colummn that ID's are held in.
  
  * {transcripts,genes}   Are you parsing through GeneIDs or Transcript IDs?

### Optional Arguments:

  * -h, --help:                 Show this help message and exit
  
  * -s START, --start START:    Row you would like to start finding variants from.
  
  * -e END, --end END:          Row you would like to stop finding variants from.

### Examples:

  ```
  python2 Variant_Finder.py Transcripts.csv ex_Transcripts_ALL.csv TranscriptID transcripts
  ```
  
  ```
  python2 Variant_Finder.py Transcripts.csv ex_Transcripts_100_to_200.csv TranscriptID transcripts -s 100 -e 200
  ```
  
  ```
  python2 Variant_Finder.py Genes.csv ex_Genes_0_to_100.csv GeneID genes -e 100
  ```
  
  ```
  python2 Variant_Finder.py Genes.csv ex_Genes_500_to_END.csv GeneID genes -s 500
  ```

# Ex_Files Folder
  * A folder that holds the output files of variants from running Variant_Finder.py
  * Fututure applcation may have folders within here for different projects.
  
# Test Folder 
  * Holds any test files to test on script.
  
  
