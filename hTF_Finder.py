"""
Aaron Rothleder
Oct 18 2016
SALK Institute PBIO-E
VariantFinder.py 
This code produces the location, variant type, allele frequency, and other data
of protein domains based on a given transcript ID.
Imports CSV file and exports CSV
"""

import requests, sys, re, pandas as pd, logging
logger = logging.getLogger('VariantFinder')
format = "%(asctime)s [%(levelname)s] %(message)s"
logging.basicConfig(filename = 'logfile.log',level=logging.INFO)

# CRITICAL - 50 - failure, application must close
# ERROR - 40 - function failed
# WARNING - 30 - unexpected ocurrance
# INFO - 20 - confirmation things went expected
# DEBUG - 10 - detailed info

#import packagges used in method
#class VariantFinder(object):
#def __init__(self, transcript):
#	self.transcript = transcript
#	#logger.info("Item created: '{}'").format(self.transcript)

def get_info(transcript):

	server = "https://rest.ensembl.org"
	ext = "/lookup/id/%s?expand=1" % (transcript)

	logging.debug("have used transcript ID for restAPI")

	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
	logging.info("if URL does not run, exit")

	transcript_list = r.json()
	logging.debug("set ensembl transcript data to transcript_list")
	logging.info("transcript_list contains all data of transcript")

	Translation_dict = transcript_list['Translation'] 
	logging.info("Translation_dict is dictionary within transcript_list containing information on translation of given transcript.")

	protein = Translation_dict['id'] 
	logging.info("access protein id that transcript is translated into")

	gene = transcript_list['Parent']
	logging.info("Ensembl gene id that transcript comes from")
	#print gene
	server = "https://rest.ensembl.org"
	ext = "/lookup/id/%s?" % (gene)
	logging.debug("access gene data")
	logging.info("get ensembl data of gene id")

	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
	logging.debug("if URL does not run, exit")

	gene_list = r.json()
	#print repr(gene_list)
	logging.info("all information of gene")
	gene_start = gene_list['start']
	logging.debug("gene_start - integer representing start of gene in terms of base pairs")
 	gene_end =  gene_list['end']
 	logging.debug("gene_start - integer representing start of gene in terms of base pairs")
 	gene_description = gene_list['description']
 	logging.info("gene_description - String description of tracription factor gene codes for")
    

	#Extract Protein information using REST API
	server = "https://rest.ensembl.org"
	ext = "/overlap/translation/%s?type=PFAM;" % (protein)
	logging.debug("protein ID used for restAPI")
	#replace the example Protein id with proteinid from given Transcript
    

	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	 
	if not r.ok:
	  r.raise_for_status()
	  sys.exit()
	#if URL does not run, exit
	 
	#r.json() is a list with one item
	protein_list = r.json()
	logging.info("protein_list is list of dictionaries contains all domains of the protein")

	pronum_list = [i for i, x in enumerate(protein_list)]
	logging.info("pronum_list is list of number of dictionaries in protein_list")

	frames = []
	for item in pronum_list:
		protein_dict = protein_list[item] 
		logging.info("turns protein_list into dictionary in order to extract individual items")
		domain = protein_dict['description'] 
		logging.debug("accessed descriptiong of protein domain")
		#protein domain (PFAM)
		start = protein_dict['start'] 
		logging.debug("accessed start aa of domain")
		#amino acid start
		end = protein_dict['end'] 
		logging.debug("accessed end aa of domain")
		#amino acid end

		server = "http://exac.hms.harvard.edu/"
		ext = "/rest/transcript/variants_in_transcript/%s" % (transcript)
		logging.debug("Use ExacBrowser REST API to access variant data")
		logging.info("http://exac.hms.harvard.edu/#rest-transcript")
		logging.info("http://exac.broadinstitute.org/about")

		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		if not r.ok:
		  r.raise_for_status()
		  sys.exit()

		Exac_list = r.json()
		logging.debug("Made Exac_list")
		logging.info("Exac_list is list of all information on variants for transcript")

		var = []
		logging.debug("created an empty list var")
		for index in Exac_list:
			var.append(((index['HGVS'],index['major_consequence'],index['allele_freq']))) 
		logging.debug("adding each variation to list var")

		pro_var = []
		for item in var:
			if 'p.' in item[0]:
				pro_var.append(item)
		logging.debug("filter protein reference sequences and put into list pro_var")
		logging.info("p. reference at http://varnomen.hgvs.org/recommendations/general/")
	
		var = []
		for item in var:
			if 'stop_gained' in item or 'frameshift_variant' in item or 'splice_donor_variant' in item or 'spice_acceptor_variant' in item:
				var.append(item)
		#filters only significant(red) consequences
		red_var=[] 
		conseq = []
		allele_freq = []
		for index in range(len(var)):
			red_var.append(var[index][0]) 
			conseq.append(var[index][1])
			allele_freq.append(var[index][2])
			#add all red variants,consequence, and frequencies to individual lists

		numbers = []
		for i in red_var:
			red = int(re.search(r'\d+', i).group())
			numbers.append(red)
			#extracts number (location) of variant

		pos_in_range = [i for i, x in enumerate(numbers) if x>=start and x<=end]
		#use pos_in_range to extract the variants in range
		red_var_in_range = []
		conseq_in_range = []
		allele_freq_in_range = []
		for i in pos_in_range:
			red_var_in_range.append(red_var[i])
			conseq_in_range.append(conseq[i])
			allele_freq_in_range.append(allele_freq[i])

		#make lists for red variants in range of domains
		variants = []
		for variant in red_var_in_range:
			variants.append(variant)

		consequences = []
		for consequence in conseq_in_range:
			consequences.append(consequence)

		freqs = []
		for frequency in allele_freq_in_range:
			freqs.append(frequency)

		if len(variants) > 0:
			number = list(range(len(variants)))
			for i in number:
				result = pd.DataFrame({'TranscriptID':[transcript],
									   'Protein Domain': [domain],
									   'Domain start (aa)':[start],
									   'Domain end (aa)':[end],
									   'variant':[variants[i]],
									   'consequence':[consequences[i]],'allele freq':[freqs[i]],
									   'Transcription Factor':transcript_list['display_name'],
									   'Gene ID':transcript_list['Parent'],
									   'Gene start (bp)':[gene_start],
									   'Gene end (bp)':[gene_end],
									   'Gene Description':[gene_description]},index = [0])
				frames.append(result)
				#creates dictionary of all data wanted for each variant
				#edit here to extract data wanted
				logging.error("Must account for null entries/no result")
		else:
			result = pd.DataFrame({'TranscriptID':[transcript],
								   'Protein Domain': [domain],
								   'Domain start (aa)':[start],
								   'Domain end (aa)':[end],
								   'variant':[variants],
								   'consequence':[consequences],'allele freq':[freqs],
								   'Transcription Factor':transcript_list['display_name'],
								   'Gene ID':transcript_list['Parent'],
								   'Gene start (bp)':[gene_start],
								   'Gene end (bp)':[gene_end],
								   'Gene Description':[gene_description]},index = [0])
			frames.append(result)
			#if there is no variant will still create empty dictionary
			#prevents concatonation error
	transcriptresult = pd.concat(frames)
	return transcriptresult	

#I/O
df = pd.read_csv('annotate_hTF_2016-11-03-all_annot.csv')
#Imports csv
print df
transcripts = []
for index, row in df.iterrows():
	transcripts.append(row['ensembl_gene_id'])
	#allows to run through list of transcripts

frames = []
for item in transcripts:
#	item_01 = VariantFinder(item)
	frames.append(get_info(item))
	#frames holds return values from method
finalresult = pd.concat(frames)
finalresult.to_csv('hTF_Final_Results.csv', index = False) 
#exports CSV