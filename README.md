This repo has everything you need to get some data on the variants within genes.

# Variant_Finder.py Running The Program (in terminal):

python2 Variant_Finder.py -h --help:


usage: Variant_Finder.py [-h] [-s START] [-e END]
                         in_File ex_File col_name {transcripts,genes}


### Positional Arguments:


  in_File               File containing ID's you would like to find variants for.
                                               
  ex_File               File exported from program with all the significant variants for provided IDs.
                                           
  col_name              Name of colummn that ID's are held in.
  
  {transcripts,genes}   Are you parsing through GeneIDs or Transcript IDs?

### Optional Arguments:

  -h, --help:                 Show this help message and exit
  
  -s START, --start START:    Row you would like to start finding variants from.
  
  -e END, --end END:          Row you would like to stop finding variants from.

### Examples:

  ```
  python2 Variant_Finder.py Transcripts.csv ex_Transcripts_ALL.csv TranscriptID transcripts
  ```
  
  ```
  python2 Variant_Finder.py Transcripts.csv ex_Transcripts_ALL.csv TranscriptID transcripts -s 100 -e 200
  ```
  
  ```
  python2 Variant_Finder.py Genes.csv ex_Genes_0_to_100.csv GeneID genes -e 100
  ```
