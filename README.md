This repo has everything you need to get some data on the variants within genes.

RUNNING THE PROGRAM (in terminal):

python2 Variant_Finder.py -h --help:

usage: Variant_Finder.py [-h] [-s START] [-e END]
                         in_File ex_File col_name {transcripts,genes}

positional arguments:
  in_File               File containing ID's you would like to find variants
                        for.
                        
  ex_File               File exported from program with all the significant
                        variants for provided IDs.
                        
  col_name              Name of colummn that ID's are held in.
  
  {transcripts,genes}   Are you parsing through GeneIDs or Transcript IDs?

optional arguments:

  -h, --help            show this help message and exit
  
  -s START, --start START
                        Row you would like to start finding variants from.
  
  -e END, --end END     Row you would like to stop finding variants from.

examples:

  python2 Variant_Finder.py Transcripts.csv ex_Transcripts_ALL.csv Transcript transcripts

  python2 Variant_Finder.py Transcripts.csv ex_Transcripts_ALL.csv Transcript transcripts -s 100 -e 200
