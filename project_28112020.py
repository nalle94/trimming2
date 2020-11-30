#!/usr/bin/env python3

import sys, argparse, gzip


###################### FUNCTIONS #########################

def check_pos(value):
	'''
	Validates if input argument is a positive integer.
	
		Parameters:
			value: The value of the input argument
		Returns:
			int_value (int): Integer of the input argument
		Raises:
			argparse.ArgumentError: If input value is not a positive integer
	'''
	#O(1)
	try:
		int_value = int(value)
		if int_value < 0 or int_value is None:
			raise argparse.ArgumentTypeError("%s is not a positive integer" % value)
	except argparse.ArgumentError as error:
		print(str(error))
		sys.exit()
	return int_value


def detect_phred(filename):
	'''
	Detects phred score for sequences in FASTQ file.
	Looks for ASCI characters in dict of the phred+33 and phred+64 encoding.
	
		Parameters:
			filename: Input FASTQ file
		Returns:
			phred_scale: Returns either 'phred+33' or 'phred+64' according to input file
		Raises:
			ValueError: If phred scale cannot be determined from input FASTQ file
	'''
	#O(n*m), where n is number of reads in the file and m is number of nucleotides per read
	#Initialize
	phred_scale = None 
	quality_data = None
	line_count = 0
	#Iterate over file until phred scale is detected
	while phred_scale == None:
		#O(n)
		for line in filename:
			line = line.strip()
			line_count += 1	
			#Identify quality score line (every 4th)
			if line_count % 4 == 0:
				#O(m)
				for char in line:
					#Look for char as keys in dict of phred+33 and phred+64 to detect phred scale
					if char in phred33 and char not in phred64:
						phred_scale = 'phred+33'
					elif char in phred64 and char not in phred33:
						phred_scale = 'phred+64'
	#If phred scale is not dectected - raise error and exit program
	if phred_scale == None:
		raise ValueError('Cannot determine phred scale from infile')
		sys.exit(1)
	return phred_scale	
				  

def trim_fixed(seq_qual, fixed_trim_5, fixed_trim_3):
	'''
	Trim fixed number of nucleotides from 5' and 3' end of reads.

		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			fixed_trim_5: Number of nucleotides to trim from 5' end of read 
			fixed_trim_3: Number of nucleotides to trim from 3' end of read
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
	'''
	#O(2m), where m is number of nucleotides per read
	(seq_trim_3, asci_trim_3, seq_trim_5, asci_trim_5) = ('', '', '', '')
	
	#Trimming from 5' end of read
	count = 0
	if fixed_trim_5:
		#O(m)
		for nuc, asci in seq_qual:
			#Trim specified number of nucleotides - and save remaining nucleotides and quality data to string variables
			if count < fixed_trim_5:
				count += 1
			else:
				seq_trim_5 += nuc
				asci_trim_5 += asci
		seq_qual = list(zip(seq_trim_5, asci_trim_5))
	
	#Trimming from 3' end of read
	count = 0
	if fixed_trim_3:
		#O(m)
		for nuc, asci in seq_qual[::-1]:
			#Trim specified number of nucleotides - and save remaining read sequence and quality data to string variables
			if count < fixed_trim_3:
				count += 1
			else:
				seq_trim_3 += nuc
				asci_trim_3 += asci

	seq_qual = list(zip(seq_trim_3[::-1], asci_trim_3[::-1]))
	return seq_qual


def trim_single_nuc_5(seq_qual, threshold):
	'''
	Trim from 5' end of read based on minimum quality of single nucleotides.
	The function trims according to specified threshold until a nucleotide with a quality above the threshold is encountered.
	
		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			threshold: Minimum phred scale quality used for trimming
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
		Raises:
			KeyError: If phred scale is incorrectly specified
	'''
	#O(m), where m is number of nucleotides per read
	(seq_trim, asci_trim) = ('', '')

	for nuc,asci in seq_qual:
		try:
			#If value of asci is below trimming threshold - trim from file until a value above threshold is encountered
			if phred[asci] < threshold and seq_trim == '':
				continue
			#Save remaining read sequence and quality data to string variables 
			else:
				seq_trim += nuc
				asci_trim += asci
		#If asci not present as key in phred scale dict - exit program
		except KeyError as error:
			sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
			sys.exit(1)	

	seq_qual = list(zip(seq_trim, asci_trim))
	return seq_qual


def trim_single_nuc_3(seq_qual, threshold):
	'''
	Trim from 3' end of read based on minimum quality of single nucleotides.
	The function trims according to specified threshold until a nucleotide with a quality above the threshold is encountered.
	
		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			threshold: Minimum phred scale quality used for trimming
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
		Raises:
			KeyError: If phred scale is incorrectly specified
	'''
	#O(m), where m is number of nucleotides per read
	(seq_trim, asci_trim) = ('', '')

	for nuc,asci in seq_qual[::-1]:
		#If value of asci is below trimming threshold - trim from file until a value above threshold is encountered
		try:
			if phred[asci] < threshold and seq_trim == '':
				continue
			#Save remaining read sequence and quality data to string variables
			else:
				seq_trim += nuc
				asci_trim += asci
		#If asci not present as key in phred scale dict - exit program
		except KeyError as error:
			sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
			sys.exit(1)	

	seq_qual = list(zip(seq_trim[::-1], asci_trim[::-1]))
	return seq_qual


def trim_moving_window_5(seq_qual, threshold):
	'''
	Trim from 5' end based on average of moving window.
	The function uses a window of size 5 and trims according to specified threshold until the average of the moving window exceeds the quality threshold.
	
		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			threshold: Minimum phred scale quality used for trimming
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
		Raises:
			KeyError: If phred scale is incorrectly specified
	'''
	#O(m), where m is number of nucleotides per read
	(seq_trim, asci_trim) = ('', '')
	window_size = 5
	i = 0

	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			sum_window = 0
			window = seq_qual[i:i + window_size]
			#If asci not present as key in phred scale dict - exit program
			try:
				for nuc,asci in window:
					sum_window += phred[asci]
			except KeyError as error:
				sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
				sys.exit(1)	
			window_average = sum_window/window_size
			#Trim 5' end of read if average of moving window is below threshold - but only until first window of satisfactory quality is encountered
			if window_average < threshold and seq_trim == '':
				i += 1
				continue
			#If average quality of moving window is ok - then save first nucleotide and quality letter to string variables and increment i
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]

		#Add remaining 3' nucleotides (last window) to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci

	#If sequence read is shorter than size of window - no trimming performed
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci
			
	seq_qual = list(zip(seq_trim, asci_trim))
	return seq_qual


def trim_moving_window_3(seq_qual, threshold):
	'''
	Trim from 3' end based on average of moving window.
	The function uses a window of size 5 and trims according to specified threshold until the average of the moving window exceeds the quality threshold.
	
		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			threshold: Minimum phred scale quality used for trimming
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
		Raises:
			KeyError: If phred scale is incorrectly specified
	'''
	#O(m), where m is number of nucleotides per read
	(seq_trim, asci_trim) = ('', '')
	window_size = 5
	i = 0

	#Read sequence and quality data from 3' end 
	seq_qual = seq_qual[::-1]
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			sum_window = 0
			window = seq_qual[i:i + window_size]
			#If asci not present as key in phred scale dict - exit program
			try:
				for nuc,asci in window:
					sum_window += phred[asci]
			except KeyError as error:
				sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
				sys.exit(1)	  
			window_average = sum_window/window_size
			#Trim 3' end of read if average of moving window is below threshold - but only until first window of satisfactory quality is encountered
			if window_average < threshold and seq_trim == '':
				i += 1
				continue
			#If average quality of moving window is ok - then save first nucleotide and quality letter to string variables and increment i
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]
			
		#Add remaining 5' nucleotides (last window) to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci

	#If sequence read is shorter than size of window - no trimming performed
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci

	seq_qual = list(zip(seq_trim[::-1], asci_trim[::-1]))
	return seq_qual


def trim_min_moving_window_5(seq_qual, threshold):
	'''Trim from 5' based on minimum of moving window.
	The function uses a window of size 5 and trims according to specified threshold until the minimum of the moving window exceeds the quality threshold.
	
		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			threshold: Minimum phred scale quality used for trimming
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
		Raises:
			KeyError: If phred scale is incorrectly specified
	'''
	#O(m), where m is number of nucleotides per read
	(seq_trim, asci_trim) = ('', '')
	window_size = 5
	i = 0

	#Read sequence and quality data from 5' end 
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			qual_window = list()
			window = seq_qual[i:i + window_size]
			#If asci not present as key in phred scale dict - exit program	
			try:
				for nuc,asci in window:
					qual_window.append(phred[asci])
			except KeyError as error:
				sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
				sys.exit(1)
			#Trim 5' end of read if minimum of moving window is below threshold - but only until first window of satisfactory quality is encountered
			if min(qual_window) < threshold and seq_trim == '':
				i += 1
				continue
			#If minimum quality of moving window is ok - then save first nucleotide and quality letter to string variables and increment i
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]

		#Add remaining 3' nucleotides (last window) to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci

	#If sequence read is shorter than size of window - no trimming performed
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci
			
	seq_qual = list(zip(seq_trim, asci_trim))
	return seq_qual


def trim_min_moving_window_3(seq_qual, threshold):
	'''Trim from 3' based on minimum of moving window.
	The function uses a window of size 5 and trims according to specified threshold until the minimum of the moving window exceeds the quality threshold.
	
		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
			threshold: Minimum phred scale quality used for trimming
		Returns:
			seq_qual: A list of tuples with paired sequence and quality data after trimming
		Raises:
			KeyError: If phred scale is incorrectly specified
	'''
	#O(m), where m is number of nucleotides per read
	(seq_trim, asci_trim) = ('', '')
	window_size = 5
	i = 0

	#Read sequence and quality data from 3' end 
	seq_qual = seq_qual[::-1]
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			qual_window = list()
			window = seq_qual[i:i + window_size]
			#If asci not present as key in phred scale dict - exit program			
			try:
				for nuc,asci in window:
					qual_window.append(phred[asci])
			except KeyError as error:
				sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
				sys.exit(1)
			#Trim 3' end of read if minimum of moving window is below threshold - but only until first window of satisfactory quality is encountered
			if min(qual_window) < threshold and seq_trim == '':
				i += 1
				continue
			#If minimum quality of moving window is ok - then save first nucleotide and quality letter to string variables and increment i
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]

		#Add remaining 5' nucleotides (last window) to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci

	#If sequence read is shorter than size of window - no trimming performed
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci
			
	seq_qual = list(zip(seq_trim[::-1], asci_trim[::-1]))
	return seq_qual


def calc_mean_qual(seq_qual):
	'''Calculate mean quality of read.

		Parameters:
			seq_qual: A list of tuples with paired sequence and quality data
		Returns:
			mean_read_qual (int): Phred scale mean quality of read
	'''
	#O(m), where m is number of nucleotides per read
	sum_qual_read = 0
	if len(seq_qual) != 0:
		for nuc, asci in seq_qual:
			#If asci not present as key in phred scale dict - exit program
			try:
				sum_qual_read += phred[asci]
			except KeyError as error:
				sys.stdout.write('Phred scale incorrectly specified, reason: ' + str(error) + ' not present in ' + str(phred_scale) + '\n')
				sys.exit(1)	
		mean_read_qual = sum_qual_read / len(seq_qual)
	#If sequence length is 0 - i.e. sequence is trimmed to length 0 - mean quality is set to 0
	else:
		mean_read_qual = 0
	return int(mean_read_qual)
	
	

############### DEFINING ARGUMENTS USING ARGPARSE ###############
	
#Read argument line	
parser = argparse.ArgumentParser(description = 'Trimming of fastq file.')
optional = parser._action_groups.pop()
required = parser.add_argument_group('Required arguments')

required.add_argument('-f', dest = 'filename', required = True, 
					help = 'input fastq filename for trimming (type: .fastq or txt.gz)')
required.add_argument('-o', dest = 'outfilename', required = True, 
					help = 'output filename. If gzip file is wanted, add .gz in end of filename (type: .fastq or txt.gz')
optional.add_argument('-p', dest = 'phred_scale', choices = ['phred+33', 'phred+64'], default = None,
					help = 'phred scale (type: phred+33 or phred+64)') 
optional.add_argument('-l', dest = 'logfile', default = 'logfile.txt', type = str,
					help = 'log filename. Need to be a txt file (default: logfile.txt)')
parser._action_groups.append(optional)


#Add mutually exclusive argument for 5' end trimming
trim_group1 = parser.add_argument_group('Trimming from 5end', description = 'Settings for trimming from 5end. Except from fixed_trim_5 - only one trimming option can be performed from each end. If no options are entered, no trimming will be performed')
trim_group1.add_argument('-t5', dest = 'fixed_trim5', default = False, type = check_pos, 
					help = 'fixed number of nucelotides to trim from 5end (type: pos int, default: False)')
trim5 = trim_group1.add_mutually_exclusive_group()
trim5.add_argument('-m5', dest = 'min_residue5', type = check_pos, default = False, 
					help = 'minimum of single residue trimming from 5end (type: pos int, default: False)') 
trim5.add_argument('-a5', dest = 'mean_mw5', type = check_pos, default = False, 
					help = 'mean of moving window trimming from 5end (type: pos int, default: False)') 
trim5.add_argument('-w5', dest = 'min_mw5', type = check_pos, default = False,
					help = 'minimum of moving window trimming from 5end (type: pos int, default: False)')


#Add mutually exclusive argument for 3' end trimming
trim_group2 = parser.add_argument_group('Trimming from 3end', description = 'Settings for trimming from 3end. Except from fixed_trim_3 - only one trimming option can be performed from each end. If no options are entered, no trimming will be performed')
trim_group2.add_argument('-t3', dest = 'fixed_trim3', default = False, type = check_pos, 
					help = 'fixed number of nucleotides to trim from 3end (type: pos int, default: False)')
trim3 = trim_group2.add_mutually_exclusive_group()
trim3.add_argument('-m3', dest = 'min_residue3', type = check_pos, default = False, 
					help = 'minimum of single residue trimming from 3end (type: pos int, default: False)') 
trim3.add_argument('-a3', dest = 'mean_mw3', type = check_pos, default = False, 
					help = 'mean of moving window trimming from 3end (type: pos int, default: False)')   
trim3.add_argument('-w3', dest = 'min_mw3', type = check_pos, default = False,
					help = 'minimum of moving window trimming from 3end (type: pos int, default: False)')


#Add grouped argument for filtering settings
filtering = parser.add_argument_group('Filtering', description = 'Settings for filtering after trimming of read')
filtering.add_argument('-q', dest = 'min_mean_qual', type = check_pos, default = 20,
					help = 'minimum mean quality after trimming (default = 20)')
filtering.add_argument('-r', dest = 'min_read_len', type = check_pos, default = 50,
					help = 'minimum read length after trimming (default = 50)')
filtering.add_argument('-s', dest = 'max_N', type = check_pos, default = 5,
					help = 'maximum number of unknown nucleotides after trimming (default = 5)')

#Save argparse arguments to a name
args = parser.parse_args()
options = vars(args)


########## READ FILES ##########

#Read input filename according to file type
try:
	if args.filename.endswith('txt.gz'):
		infile = gzip.open(args.filename,'rt')
	elif args.filename.endswith('.fastq'):
		infile = open(args.filename, 'r')
except IOError as error:
		sys.stdout.write('Cannot open inputfile, reason: ' + str(error) + '\n')
		sys.exit(1)

#Read output filename according to file type
try:
	if args.outfilename.endswith('txt.gz'):
		outfile = gzip.open(args.outfilename, 'wt')
	else:
		outfile = open(args.outfilename, 'w')
except IOError as error:
	sys.stdout.write('Cannot open outfile, reason: ' + str(error) + '\n')
	sys.exit(1)
  

#Read logfile
if args.logfile.endswith('.txt'):
	logfile = open(args.logfile, 'w')
else:
	logfile = open('logfile.txt', 'w')

#Add header to logfile
print('{:15}'.format('Original file:'), args.filename, '\n''{:15}'.format('Trimmed file:'), args.outfilename, '\n', file = logfile)



########### DIFFERENT PHRED SCORE ENCODINGS ###########

phred33 = {
	'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, 
	'+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19,
	'5': 20, '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29,
	'?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39,
	'I': 40, 'J': 41}
phred64 = {
	'@': 0, 'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8,	'I': 9,
	'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 'O': 15, 'P': 16, 'Q': 17, 'R': 18, 'S': 19,
	'T': 20, 'U': 21, 'V': 22, 'W': 23, 'X': 24, 'Y': 25, 'Z': 26, '[': 27, '\\': 28, ']': 29,
	'^': 30, '_': 31, '`': 32, 'a': 33, 'b': 34, 'c': 35, 'd': 36, 'e': 37, 'f': 38, 'g': 39, 
	'h': 40, 'i': 41}


########## MAIN PROGRAM ##########

#Detect phred score encoding from input file - if not given as input argument
if args.phred_scale == None:
	#O(n*m), where n is number of reads in the file and m is number of nucleotides per read
	phred_scale = detect_phred(infile)
	#Reset position in file to beginning of file
	infile.seek(0)
else:
	phred_scale = args.phred_scale
	
#Set phred scale encoding
if phred_scale == 'phred+33':
	phred = phred33
elif phred_scale == 'phred+64':
	phred = phred64


#Initialize
(read_count, trim_count, filtered_count, line_count, nuc_sum) = (0, 0, 0, 0, 0)
(header, seq, plus, qual, seq_qual_trim) = ('', '', '', '', '')
(count_n, count_a, count_c, count_g, count_t, count_n_trim)= (0, 0, 0, 0, 0, 0)
seq_length = list()

#Read one sequence at a time - 4 lines
#O(n), where n is number of reads in the file
for line in infile:
	line = line.strip()
	line_count += 1
	#Identify header, sequence, third line and quality data for read
	try:
		if line_count % 4 == 1:
			header = line
			#If header does not start with @ make error in filetype
			if header.startswith('@') is not True:
				raise IOError('Incorrect header in input file read, ' + header + ', header is required to start with "@"')
		elif line_count % 4 == 2:
			seq = line
		elif line_count % 4 == 3:
			plus = line
			#If 3. line does not start with "+"" make error in filetype
			if plus.startswith('+') is not True:
				raise IOError('Incorrect third line in input file read, ' + header + ', third line is required to start with "+"')
		elif line_count % 4 == 0:
			qual = line

		#If sequence and quality data are both containing data - trimming functions are performed according to arguments from command line
		if seq != '' and qual != '':
			#Set phred scale encoding
			if phred_scale == 'phred+33':
				phred = phred33
			elif phred_scale == 'phred+64':
				phred = phred64

			#Create list of tuples - where each sequence letter is combined with its phred scale quality letter
			seq_qual = list(zip(seq, qual))	
			#Add read sequence length to list
			seq_length.append(len(seq))
			#Count number of nucleotides in input file
			for nuc, asci in seq_qual:
				count_n += nuc.count('N')
				count_a += nuc.count('A')
				count_c += nuc.count('C')
				count_g += nuc.count('G')
				count_t += nuc.count('T')
			"""Check if sequence only contains A, T, C, G or N by comparing the total nucleotide count
			with the previous total nucleotide count plus the length of sequence"""
			if (len(seq) + nuc_sum) != count_a + count_t + count_c + count_g + count_n:
				raise IOError('Wrong symbol in sequence line of read, ' + header + ', only A, C, T, G and N are allowed')
			#Add length of sequence to total nucleotide count
			nuc_sum += len(seq)
			#Check if length of sequence and quality data is identical
			if len(seq) != len(qual):
				raise IOError('Length of sequence and quality data does not match in read, ' + header)

			seq_qual_trim = seq_qual[:]

			#Trim fixed number of nucleotides - O(2m), where m is number of nucleotides per read
			if args.fixed_trim5 or args.fixed_trim3:
				seq_qual_trim = trim_fixed(seq_qual_trim, args.fixed_trim5, args.fixed_trim3)
			#Trim based on quality of nucleotides from 5' end - O(m), where m is number of nucleotides per read
			if args.min_residue5:
				seq_qual_trim = trim_single_nuc_5(seq_qual_trim, args.min_residue5)
			elif args.mean_mw5:
				seq_qual_trim = trim_moving_window_5(seq_qual_trim, args.mean_mw5)
			elif args.min_mw5:
				seq_qual_trim = trim_min_moving_window_5(seq_qual_trim, args.min_mw3)
			#Trim based on quality of nucleotides from 3' end  - O(m), where m is number of nucleotides per read
			if args.min_residue3:
				seq_qual_trim = trim_single_nuc_3(seq_qual_trim, args.min_residue3)
			elif args.mean_mw3:
				seq_qual_trim = trim_moving_window_3(seq_qual_trim, args.mean_mw3)
			elif args.min_mw3:
				seq_qual_trim = trim_min_moving_window_3(seq_qual_trim, args.min_mw3)

			#Increase read count by 1
			read_count += 1

			#Count number of trimmed reads
			if len(seq_qual_trim) < len(seq_qual):
				trim_count += 1

			#Filter reads based on users input - O(m), where m is number of nucleotides per read
			read_mean_qual = calc_mean_qual(seq_qual_trim)	 		#Calculate mean quality of read after trimming 
			for nuc, asci in seq_qual_trim:							#Calculate number of unknown nucleotides in read after trimming
				count_n_trim += nuc.count('N')

			if read_mean_qual < args.min_mean_qual or len(seq_qual_trim) < args.min_read_len or count_n_trim > args.max_N:
				filtered_count += 1
				#Reset variables and continue to next read with no saving to outfile
				(header, seq, plus, qual, count_n_trim) = ('', '', '', '', 0)
				continue

			#Print read to outfile
			seq_qual_sep = list(zip(*seq_qual_trim))				#Create list with two lists - sequence and phred scale quality letters, respectively
			print(header, file = outfile)
			print(''.join(seq_qual_sep[0]), file = outfile)
			print(plus, file = outfile)
			print(''.join(seq_qual_sep[1]), file = outfile)
			
			#Reset variables
			(header, seq, plus, qual, count_n_trim) = ('', '', '', '', 0)

	except IOError as error:
		sys.stdout.write('Incorrect input file format, reason: ' + str(error) + '\n')
		sys.exit(1)

		
#Saving data to logfile
print('{:25}'.format('Settings for analysis'), file = logfile)
[print(key, ':', value, file = logfile) for key, value in options.items()]
print('\nCounts of reads in file', file = logfile)
print('{:25}'.format('Number of reads in file:'), '{:>10}'.format(read_count), file = logfile)
print('{:25}'.format('Number of trimmed reads:'), '{:>10}'.format(trim_count), file = logfile)
print('{:25}'.format('Number of filtered reads:'), '{:>10}'.format(filtered_count), file = logfile)

#Calculate and print nucleotide counts and GC content in input file
GC = (float(count_g) + float(count_c)) / float(nuc_sum)
print('\nTotal nucleotide counts in input file\nA:', count_a, '\nT:', count_t, '\nG:', count_g, '\nC:', count_c, '\nGC content:', '{:.2f}'.format(GC), file = logfile)

#Print information on read lengths in input file
print('\nAverage read length in input file: {0}'.format(sum(seq_length)/len(seq_length)), file = logfile)
print('Minimum read length in input file: {0}'.format(min(seq_length)), file = logfile)
print('Maximum read length in input file: {0}'.format(max(seq_length)), file = logfile)

#Close files
infile.close()
outfile.close()
logfile.close()



