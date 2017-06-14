"""
Aaron Rothleder
Oct 18 2016
SALK Institute PBIO-E
VariantFinder.py 
This code produces the location, variant type, allele frequency, and other data
of protein domains based on a given transcript ID.
Imports CSV file and exports CSV
"""
##########################Logging Info##########################
# CRITICAL - 50 - failure, application must close
# ERROR - 40 - function failed
# WARNING - 30 - unexpected ocurrance
# INFO - 20 - confirmation things went expected
# DEBUG - 10 - detailed info
################################################################

import requests, sys, re, pandas as pd, logging, urllib2, time, argparse

#Argument Parsing
parser = argparse.ArgumentParser()

parser.add_argument("in_File", 
					help="File containing ID's you would like to find variants for.")
parser.add_argument("ex_File",
					help="File exported from program with all the significant variants for provided IDs.")
parser.add_argument("col_name",
					help="Name of colummn that ID's are held in.")
parser.add_argument("ID_type", choices = ["transcripts", "genes"],
					help="Are you parsing through GeneIDs or Transcript IDs?")
parser.add_argument("-s", "--start", type=int,
					help="Row you would like to start finding variants from.")
parser.add_argument("-e", "--end", type=int,
					help="Row you would like to stop finding variants from.")

args = parser.parse_args()


#main method
def main():
	print str(sys.argv)

	#in_File, Ex_file, cole_name, ID_type assignment
	print "Reading file " + args.in_File
	in_File = args.in_File
	print "Will export file " + args.ex_File
	ex_File = args.ex_File
	print "IDs should be in column " + args.col_name
	col_name = args.col_name
	print "IDs should be " + args.ID_type
	ID_type = args.ID_type


	
	"""COUNTERS:"""
	geneNumber = 0
	#global transcriptNumber 
	transcriptNumber = 0
	
	#I/O
	totalStart = time.time()
	Largedf = pd.read_csv(in_File)
	#print type(start)
	#print end
	#print "PRINTING LARGEdf"
	#print Largedf
	#Smalldf = Largedf.head(5)
	if args.start and not args.end:
		start = sys.args.start
		end = len(Largedf)
	elif not args.start and args.end:
		start = 0
		end = args.end
	elif args.start and args.end:
		if args.start >= args.end:
			print "End must be greater than start"
			exit()
		start = args.start
		end = args.end
	else:
		start = 0
		end = len(Largedf)

	print "Running from {} to {} .".format(start,end)
	Smalldf = Largedf[int(start):int(end)]
	#use start and end to create chunk
	#Imports csv
	print "PRINTING what is will be ran through script."
	print Smalldf
	#print Smalldf
	genes = []
	transcripts_frames=[]

	if ID_type == "genes":
		for index, row in Smalldf.iterrows():
			#allows to run through list of transcripts
			genes.append(row[col_name])
		#print "Printing genes"
		#print genes
	


		geneNumber = len(genes)
		for item in genes:
			geneNumber += 1
			start = time.time()
			server = "http://exac.hms.harvard.edu/"
			ext = "/rest/gene/%s" % (item)
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			end = time.time()
			timePassed = end-start
			print "API: Time to get gene info through Exac: " + str(timePassed)
			
			
			start = time.time()
			try:
				gene_info = r.json()
				#if gene_info is not None:
				transcripts_in_gene = gene_info['transcripts_in_gene']
				end = time.time()
				timePassed = end-start
				print "Time to get transcripts_in_gene info: " + str(timePassed)
				transcripts = []
				#transcripts_frames = []
				for index in range(len(transcripts_in_gene)):
					transcripts.append(transcripts_in_gene[index]['transcript_id'])
				for index in range(len(transcripts)):
					transcripts_frames.append(transcripts[index])
			except KeyError as err:
				loggin.warning("Problem with accessing Exac")
			except ValueError as err:
				logging.warning("gene not found in Exac")

	elif ID_type == "transcripts":
		for index, row in Smalldf.iterrows():
			transcripts_frames.append(row[col_name])

	frames = []
	empty = []
	for item in transcripts_frames:
		print "still running"
		#print item
		try:
			#call to get_info
			check = get_info(item)
			#item_01 = VariantFinder(item)
			frames.append(check)
			#frames holds return values from method
			finalresult = pd.concat(frames)
			finalresult.to_csv(ex_File, index = False)


		except ValueError as err:
			logging.warning("item not in dictionary")
		except KeyError as err:
			logging.warning("item not in dictionary")
		except urllib2.HTTPError as err:
			loggng.warning("incorrect URL")

	print "method finished!" 
	totalEnd = time.time()
	totalTimePassed = totalEnd - totalStart
	print str(geneNumber) + " genes | " + str(transcriptNumber) + " transcripts"
	print "Total program runtime: " + str(totalTimePassed)
	#exports CSV
		
		
#get_info method
def get_info(transcript):	
	logger = logging.getLogger('VariantFinder')
	ch = logging.StreamHandler()
	ch.setLevel(logging.DEBUG)
	format = "%(asctime)s [%(levelname)s] %(message)s"
	logging.basicConfig(filename = 'logfile.log',level=logging.WARNING)
	logger.addHandler(ch)


	server = "https://rest.ensembl.org"
	ext = "/lookup/id/%s?expand=1" % (transcript)

	#transcriptNumber += 1
	print transcript

	logging.debug("have used transcript ID for restAPI")
	start = time.time()
	r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
	end = time.time()
	timePassed = end-start
	print "API: ensembl transcript lookup/id time: " + str(timePassed)


	#r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
 	"""
	if not r.ok:
	  r.raise_for_status()
	  sys.close()
	"""

	logging.info("if URL does not run, exit")


	transcript_list = r.json()
	#print transcript_list

	logging.debug("set ensembl transcript data to transcript_list")
	logging.info("transcript_list contains all data of transcript")
	biotype = transcript_list['biotype']
	#print biotype

	#####################################Protein Method#########################
	if biotype == "protein_coding":
		Translation_dict = transcript_list['Translation'] 

		logging.info("Translation_dict is dictionary within transcript_list containing information on translation of given transcript.")

		protein = Translation_dict['id'] 
		logging.info("access protein id that transcript is translated into")
		#print protein

		gene = transcript_list['Parent']
		logging.info("Ensembl gene id that transcript comes from")
		#print gene
		start = time.time()
		server = "https://rest.ensembl.org"
		ext = "/lookup/id/%s?" % (gene)
		logging.debug("access gene data")
		logging.info("get ensembl data of gene id")

		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		end = time.time()
		timePassed = end - start
		print "API: ensembl gene data load time: " + str(timePassed)

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

	 	#################################Variants in Transcript Method#####################
		start = time.time()
		server = "http://exac.hms.harvard.edu/"
		ext = "/rest/transcript/variants_in_transcript/%s" % (transcript)
		logging.debug("Use ExacBrowser REST API to access variant data")
		logging.debug("http://exac.hms.harvard.edu/#rest-transcript")
		logging.debug("http://exac.broadinstitute.org/about")

		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
		end = time.time()
		timePassed = end-start
		print "API: time to get transcipt info from Exac: " + str(timePassed)

		if not r.ok:
		  r.raise_for_status()
		  sys.exit()

		Exac_list = r.json()
		#print Exac_list
		logging.debug("Made Exac_list")
		logging.info("Exac_list is list of all information on variants for transcript")

		var = []
		logging.debug("created an empty list var")
		for index in Exac_list:
			var.append(((index['HGVS'],index['major_consequence'],index['allele_freq']))) 
		logging.debug("adding each variation to list var")

		pro_var = []
		for item in var:
			if 'p.' or 'c.' in item[0]:
				pro_var.append(item)
		logging.debug("filter protein reference sequences and put into list pro_var")
		logging.info("p. reference at http://varnomen.hgvs.org/recommendations/general/")
	
		var = []
		for item in pro_var:
			if 'stop_gained' in item or 'frameshift_variant' in item or 'splice_donor_variant' in item or 'splice_acceptor_variant' in item or 'synonymous_variant' in item:
				
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
	    

########################################Protein Method#################################
		#Extract Protein information using REST API
		start = time.time()
		server = "https://rest.ensembl.org"
		ext = "/overlap/translation/%s?;" % (protein)
		logging.debug("protein ID used for restAPI")
		#replace the example Protein id with proteinid from given Transcript
	    

		r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

		end = time.time()
		timePassed = end - start
		print "API: load time of protein info using ensembl: " + str(timePassed)
		
		#if URL does not run, exit
		 
		#r.json() is a list with one item
		protein_list = r.json()

		logging.info("protein_list is list of dictionaries contains all domains of the protein")

		pronum_list = [i for i, x in enumerate(protein_list)]
		logging.info("pronum_list is list of number of dictionaries in protein_list")

		if len(pronum_list)>0:
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

#########################################Output Method################################
				if len(variants) > 0:
					number = list(range(len(variants)))
					for i in number:
						result = pd.DataFrame({'TranscriptID':[transcript],
											   'Protein Domain': [domain],
											   'Domain start (aa)':[start],
											   'Domain end (aa)':[end],
											   'variant':[variants[i]],
											   'consequence':[consequences[i]],
											   'allele freq':[freqs[i]],
											   'Transcription Factor':transcript_list['display_name'],
											   'Gene ID':transcript_list['Parent'],
											   'Gene start (bp)':[gene_start],
											   'Gene end (bp)':[gene_end],
											   'Gene Description':[gene_description],
											   'Domain length (aa)':[end-start]},index = [0])
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
										   'consequence':[consequences],
										   'allele freq':[freqs],
										   'Transcription Factor':transcript_list['display_name'],
										   'Gene ID':transcript_list['Parent'],
										   'Gene start (bp)':[gene_start],
										   'Gene end (bp)':[gene_end],
										   'Gene Description':[gene_description],
										   'Domain length (aa)':[end-start]},index = [0])
					frames.append(result)
					#if there is no variant will still create empty dictionary
					#prevents concatonation error
			transcriptresult = pd.concat(frames)
			return transcriptresult




#main method call
if __name__ == "__main__":
	main()