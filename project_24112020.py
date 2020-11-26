#!/usr/bin/env python3

import sys, argparse, gzip


###################### FUNCTIONS #########################

def check_pos(value):
	'''validate if args is a pos int'''
	try:
		int_value = int(value)
		if int_value < 0:
			raise argparse.ArgumentError("%s is not a positive int value" % value)
	except argparse.ArgumentError as error:
		print("%s is not a positive int value" % value)
		sys.exit()
	return int_value

def check_qual_cut(value):
	'''check if value is pos int or return defaultvalue if yes'''
	if value == 'yes':
		value = 0
	else:
		try:
			int_value = int(value)
			if int_value < 0:
				raise argparse.ArgumentError("%s is not a positive int value" % value)
		except argparse.ArgumentError as error:
			print("%s is not a positive int value" % value)
			sys.exit()
	return int_value


def detect_phred(filename):
	'''detects phred score for sequences in file'''
	#Set start
	phred_scale = None 
	quality_data = None
	line_count = 0
	#iterate over file until scale is identified
	while phred_scale == None:
		for line in filename:
			line = line.strip()
			line_count += 1	
			#identify quality score lines (every 4th)
			if line_count % 4 == 0:
				for char in line:
					#translate the line to find phred scales
					if char in phred33 and char not in phred64:
						phred_scale = "phred+33"
					elif char in phred64 and char not in phred33:
						phred_scale = "phred+64"
	#if program can not determine phred scale exit program
	if phred_scale == None:
		print("Cannot determine phred scale")
		sys.exit(1)
	return phred_scale	
				  

def trim_fixed(seq_qual, fixed_trim_5, fixed_trim_3):
	'''Trim fixed number of nucleotides from 5' and 3' ends'''
	count = 0
	seq_trim_3 = ''
	asci_trim_3 = ''
	seq_trim_5 = ''
	asci_trim_5 = ''
	#if trimming from 5end is chosen
	if fixed_trim3:
		for nuc, asci in seq_qual:
			if count < fixed_trim_5:
				count += 1
			else:
				seq_trim_5 += nuc
				asci_trim_5 += asci
	seq_qual_trim_5 = list(zip(seq_trim_5, asci_trim_5))
	count = 0
	#if trimming from 3end is chosen
	if fixed_trim3:
		for nuc, asci in seq_qual_trim_5[::-1]:
			if count < fixed_trim_3:
				count += 1
			else:
				seq_trim_3 += nuc
				asci_trim_3 += asci

	seq_qual = list(zip(seq_trim_3[::-1], asci_trim_3[::-1]))
	return seq_qual


def trim_single_nuc_5(seq_qual, threshold):
	'''Trim from 5' based on minimum quality of single nucleotides'''
	seq_trim_5 = ''
	asci_trim_5 = ''
	for nuc,asci in seq_qual:
		#if asci not present in phred scale exit
		try:
			if phred[asci] < threshold and seq_trim_5 == '':
				continue
			else:
				seq_trim_5 += nuc
				asci_trim_5 += asci
		except KeyError as error:
			sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
			sys.exit(1)	

	seq_qual = list(zip(seq_trim_5, asci_trim_5))
	return seq_qual


def trim_single_nuc_3(seq_qual, threshold):
	'''Trim from 3' based on minimum quality of single nucleotides'''
	seq_trim_3 = ''
	asci_trim_3 = ''
	for nuc,asci in seq_qual[::-1]:
		#if asci not present in phred scale exit
		try:
			if phred[asci] < threshold and seq_trim_3 == '':
				continue
			else:
				seq_trim_3 += nuc
				asci_trim_3 += asci
		except KeyError as error:
			sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
			sys.exit(1)	
	seq_qual = list(zip(seq_trim_3[::-1], asci_trim_3[::-1]))
	return seq_qual


def trim_moving_window_5(seq_qual, threshold):
	'''Trim from 5' based on average of moving window of size 5'''
	seq_trim = ''
	asci_trim = ''
	window_size = 5
	i = 0
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			sum_window = 0
			window = seq_qual[i:i + window_size]
			for nuc,asci in window:
				#if asci not present in phred scale, exit
				try:
					sum_window += phred[asci]
				except KeyError as error:
					sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
					sys.exit(1)	
			window_average = sum_window/window_size
			if window_average < threshold and seq_trim == '':
				i += 1
				continue
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]			
		#Add remaining 3' nucleotides to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci
	#if there is less nucleotides left than the size of window
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci
			
	seq_qual = list(zip(seq_trim, asci_trim))
	return seq_qual


def trim_moving_window_3(seq_qual, threshold):
	'''Trim from 3' based on average of moving window of size 5'''
	seq_trim = ''
	asci_trim = ''
	window_size = 5
	i = 0
	seq_qual = seq_qual[::-1]
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			sum_window = 0
			window = seq_qual[i:i + window_size]
			#if asci not present in phred scale, exit
			try:
				for nuc,asci in window:
					sum_window += phred[asci]
			except KeyError as error:
				sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
				sys.exit(1)	  
			window_average = sum_window/window_size
			if window_average < threshold and seq_trim == '':
				i += 1
				continue
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]
			
		#Add remaining 5' nucleotides to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci
	#if there is less nucleotides left than the size of window
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci

	seq_qual = list(zip(seq_trim[::-1], asci_trim[::-1]))
	return seq_qual


def trim_min_moving_window_5(seq_qual, threshold):
	'''Trim from 5' based on minimum of moving window of size 5'''
	seq_trim = ''
	asci_trim = ''
	window_size = 5
	i = 0
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			qual_window = list()
			window = seq_qual[i:i + window_size]
			for nuc,asci in window:
				#if asci not present in phred scale, exit
				try:
					qual_window.append(phred[asci])
				except KeyError as error:
					sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
					sys.exit(1)	
			if min(qual_window) < threshold and seq_trim == '':
				i += 1
				continue
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]
		#Add remaining 3' nucleotides to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci
	#if there is less nucleotides left than the size of window
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci
			
	seq_qual = list(zip(seq_trim, asci_trim))
	return seq_qual


def trim_min_moving_window_3(seq_qual, threshold):
	'''Trim from 3' based on minimum of moving window of size 5'''
	seq_trim = ''
	asci_trim = ''
	window_size = 5
	i = 0
	seq_qual = seq_qual[::-1]
	if len(seq_qual) >= window_size:
		while i < len(seq_qual) - window_size + 1:
			qual_window = list()
			window = seq_qual[i:i + window_size]
			for nuc,asci in window:
			#if asci not present in phred scale, exit
				try:
					qual_window.append(phred[asci])
				except KeyError as error:
					sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
					sys.exit(1)	
			if min(qual_window) < threshold and seq_trim == '':
				i += 1
				continue
			else:
				i += 1
				seq_trim += window[0][0]
				asci_trim += window[0][1]

		#Add remaining 5' nucleotides to trimmed sequence
		for nuc,asci in window[1:]:
			seq_trim += nuc
			asci_trim += asci
	#if there is less nucleotides left than the size of window
	else:
		for nuc,asci in seq_qual:
			seq_trim += nuc
			asci_trim += asci
			
	seq_qual = list(zip(seq_trim, asci_trim))
	return seq_qual


def calc_mean_qual(seq_qual):
	'''Calculate mean quality of read'''
	sum_qual_read = 0
	if len(seq_qual) != 0:
		for nuc, asci in seq_qual:
			#if asci not present in phred scale exit
			try:
				sum_qual_read += phred[asci]
			except KeyError as error:
				sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
				sys.exit(1)	
		average_read_qual = sum_qual_read / len(seq_qual)
	#if there is no qual data the mean is 0
	else:
		average_read_qual = 0
	return int(average_read_qual)
	
	

###################### ARGUMENT LINE ######################
	
#read argument line	
parser = argparse.ArgumentParser(description = 'Trimming of fastq file.')
parser.add_argument('-f', dest = 'filename', required = True, 
					help = 'input fastq filename for trimming (type: fastq or fastq.gz)')
parser.add_argument('-o', dest = 'outfilename', required = True, 
					help = 'output filename. If gzip file is wanted, add .gz in end of filename (type: fastq or fastq.gz')
parser.add_argument('-p', dest = 'phred_scale', choices = ['phred+33', 'phred+64'], default = None,
					help = 'phred scale (type: phred+33 or phred+64)') 
parser.add_argument('-l', dest = 'logfile', default = 'logfile.txt', type = str,
					help = 'log filename. Need to be a txt file (default: logfile.txt)')


#add mutually exclusive argument for 3' end trimming
trim_group1 = parser.add_argument_group('Trimming from 3end', description = 'Settings for trimming from 3end. Only one trimming option can be performed from each end. If no options are entered, no trimming will be performed')
trim_group1.add_argument('-t3', dest = 'fixed_trim3', default = False, type = check_pos, 
					help = 'what fixed base length to trim (type: pos int)')
trim3 = trim_group1.add_mutually_exclusive_group()
trim3.add_argument('-m3', dest = 'min_residue3', type = check_pos, default = False, 
					help = 'minimum of single residue trimming from 3end (type: pos int, default: False)') 
trim3.add_argument('-a3', dest = 'mean_mw3', type = check_pos, default = False, 
					help = 'mean of moving window trimming from 3end (type: pos int, default: False)')   
trim3.add_argument('-w3', dest = 'min_mw3', type = check_pos, default = False,
					help = 'minimum of moving window trimming from 3end (type: pos int, default: False)')



#add mutually exclusive argument for 5' end trimming
trim_group2 = parser.add_argument_group('Trimming from 5end', description = 'Settings for trimming from 5end. Only one trimming option can be performed from each end. If no options are entered, no trimming will be performed')
trim_group2.add_argument('-t5', dest = 'fixed_trim5', default = False, type = check_pos, 
					help = 'what fixed base length to trim (type: pos int)')
trim5 = trim_group2.add_mutually_exclusive_group()
trim5.add_argument('-m5', dest = 'min_residue5', type = check_pos, default = False, 
					help = 'minimum of single residue trimming from 5end (type: pos int, default: False)') 
trim5.add_argument('-a5', dest = 'mean_mw5', type = check_pos, default = False, 
					help = 'mean of moving window trimming from 5end (type: pos int, default: False)') 
trim5.add_argument('-w5', dest = 'min_mw5', type = check_pos, default = False,
					help = 'minimum of moving window trimming from 5end (type: pos int, default: False)')

#add grouped argument for filtering settings
filtering = parser.add_argument_group('Filtering', description = 'Settings for filtering (default: min_mean_qual = 20, min_read_len = 50, max_N = 5)')
filtering.add_argument('-q', dest = 'min_mean_qual', type = check_pos, default = 20,
					help = 'minimum mean qual after trimming (default = 20)')
filtering.add_argument('-r', dest = 'min_read_len', type = check_pos, default = 50,
					help = 'minimum read length after trimming (default = 50)')
filtering.add_argument('-s', dest = 'max_N', type = check_pos, default = 5,
					help = 'maximum unknown N after trimming (default = 5)')

#save argparse arguments to a name
args = parser.parse_args()
options = vars(args)


###################### READ FILES ###########################

#read input file if filetype is gzipped fastq 
if args.filename.endswith('txt.gz'):
	try:
		infile = gzip.open(args.filename,'rt')
	except IOError as error:
		sys.stdout.write('Cannot open inputfile, reason: ' + str(error) + '\n')
		sys.exit(1)
#read file if filetype is fastq
elif args.filename.endswith('.fastq'):
	try:
		infile = open(args.filename, 'r')
	except IOError as error:
		sys.stdout.write('Cannot open inputfile, reason: ' + str(error) + '\n')
		sys.exit(1)
else:
	print('Not a valid filename')

#open output file and check if name if gzip
if args.outfilename.endswith('txt.gz'):
	try:
		outfile = gzip.open(args.outfilename, 'wt')
	except IOError as error:
		sys.stdout.write('Cannot open outfile, reason: ' + str(error) + '\n')
		sys.exit(1)
else:
	try:
		outfile = open(args.outfilename, 'w')
	except IOError as error:
		sys.stdout.write('Cannot open outfile, reason: ' + str(error) + '\n')
		sys.exit(1)   

#read logfile
if args.logfile.endswith('.txt'):
	try:
		logfile = open(args.logfile, 'w')
	except IOError as error:
		sys.stdout.write('Cannot open logfile, reason: ', + str(error) + '\n')
		sys.exit()
#add text to logfile
print('{:15}'.format('Original file:'), args.filename, '\n''{:15}'.format('Trimmed file:'), args.outfilename, '\n', file = logfile)



#### 2. Detect Phred score ####
phred33 = {
	'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6, '(': 7, ')': 8, '*': 9, 
	'+': 10, ',': 11, '-': 12, '.': 13, '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19,
	'5': 20, '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27, '=': 28, '>': 29,
	'?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34, 'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39,
	'I': 40, 'J': 41}
phred64 = {
	'@': 0, 'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8,	'I': 9,
	'J': 10, 'K': 11, 'L': 12, 'M': 13, 'N': 14, 'O': 15, 'P': 16, 'Q': 17, 'R': 18, 'S': 19,
	'T': 20, 'U': 21, 'V': 22, 'W': 23, 'X': 24, 'Y': 35, 'Z': 26, '[': 27, '\\': 28, ']': 29,
	'^': 30, '_': 31, '`': 32, 'a': 33, 'b': 34, 'c': 35, 'd': 36, 'e': 37, 'f': 38, 'g': 39, 
	'h': 40, 'i': 41}


	

####MAIN PROGRAM####


#if no phred scale is given, detect it automaticly
if args.phred_scale == None:
	phred_scale = detect_phred(infile)
	#go back to beginning of file
	infile.seek(0)
 #if phred scale is given
else:
	phred_scale = args.phred_scale
	
 
#Initialize
(read_count, trim_count, filtered_count, excluded_count, line_count, nuc_sum) = (0, 0, 0, 0, 0,0)
(header, seq, plus, qual, seq_qual_trim) = ('', '', '', '','')
(count_n, count_a, count_c, count_g, count_t, count_n_trim)= (0, 0, 0, 0, 0, 0)
error_seq = []

#read one sequence at a time and sort lines
for line in infile:
	line = line.strip()
	line_count += 1
	#initialize values
	#Identify header, sequence and quality data for read
	if line_count % 4 == 1:
		header = line
		#if header does not start with @ make error in filetype
		try:
			True == header.startswith('@') 
		except IOerror as error:
			sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
			sys.exit(1)
	elif line_count % 4 == 2:
		seq = line
	elif line_count % 4 == 3:
		plus = line
		#if 3. line does not start with + make error in filetype
		try:
			True == plus.startswith('+') 
		except IOerror as error:
			sys.stdout.write('Incorrect inputfile, reason: ' + str(error) + '\n')
			sys.exit(1)
	elif line_count % 4 == 0:
		qual = line

	#If sequence and quality data are both containing data - trimming functions are performed according to arguments from command line
	if seq != '' and qual != '':
		#naem phred scale
		if phred_scale == 'phred+33':
			phred = phred33
		elif phred_scale == 'phred+64':
			phred = phred64
		seq_qual = list(zip(seq, qual))	
		#Count number of nucleotides in input file
		for nuc, asci in seq_qual:
			count_n += nuc.count('N')
			count_a += nuc.count('A')
			count_c += nuc.count('C')
			count_g += nuc.count('G')
			count_t += nuc.count('T')
		#error check if sequence only contains A, T, C, G or N by calculating the sum and see if it matchen len of line
		if (len(seq) + nuc_sum) != count_a + count_t + count_c + count_g + count_n:
			raise IOError('Wrong symbol in sequnence line of read: ' + header + ', only A, C, T, G and N are allowed\n')
			sys.exit(1)
		#make new old nucleotide count
		nuc_sum += len(seq)
		#If length of seq and qual data is not the same, the read is not trimmed or saved to outfile
		if len(seq) != len(qual):
			raise IOError('Length of sequence and quality data does not match in read: ' + header)
	
		#Trim fixed number of nucleotides
		if args.fixed_trim5 or args.fixed_trim3:
			seq_qual_trim = trim_fixed(seq_qual, args.fixed_trim5, args.fixed_trim3)
		#Trim based on quality of nucleotides from 5' end
		if args.min_residue5:
			seq_qual_trim = trim_single_nuc_5(seq_qual_trim, args.min_residue5)
		elif args.mean_mw5:
			seq_qual_trim = trim_moving_window_5(seq_qual_trim, args.mean_mw5)
		elif args.min_mw5:
			seq_qual_trim = trim_min_moving_window_5(seq_qual_trim, args.min_mw3)
		#Trim based on quality of nucleotides from 3' end
		if args.min_residue3:
			seq_qual_trim = trim_single_nuc_3(seq_qual_trim, args.min_residue3)
		elif args.mean_mw3:
			seq_qual_trim = trim_moving_window_3(seq_qual_trim, args.mean_mw3)
		elif args.min_mw3:
			seq_qual_trim = trim_min_moving_window_3(seq_qual_trim, args.min_mw3)

		#Increase read count by 1
		read_count += 1
		#Count number of trimmed reads
		if seq_qual_trim == '':
			trim_count += 1

		#Filter reads based on users input
		read_mean_qual = calc_mean_qual(seq_qual_trim)	 	#calculate mean qual of read after trimming 
		for nuc, asci in seq_qual_trim:								#calculate number of unknown nucleotides in read after trimming
			count_n_trim += nuc.count('N')

		if read_mean_qual < args.min_mean_qual or len(seq_qual_trim) < args.min_read_len or count_n_trim > args.max_N:
			filtered_count += 1
			#reset abd continue to next read with no saving of read to outfile
			(header, seq, plus, qual, seq_qual_trim) = ('', '', '', '','')
			continue

		#Print read to outfile
		seq_qual_sep = list(zip(*seq_qual_trim))
		print(header, file = outfile)
		print(''.join(seq_qual_sep[0]), file = outfile)
		print(plus, file = outfile)
		print(''.join(seq_qual_sep[1]), file = outfile)
		
		#reset 
		(header, seq, plus, qual, seq_qual_trim) = ('', '', '', '','')

####Saving to log file
print('{:25}'.format('Settings for analysis'), file = logfile)
[print(key, ':', value, file = logfile) for key, value in options.items()]
print('\nCounts of reads in file', file = logfile)
print('{:25}'.format('Number of reads in file:'), '{:>10}'.format(read_count), file = logfile)
print('{:25}'.format('Number of trimmed reads:'), '{:>10}'.format(trim_count), file = logfile)
print('{:25}'.format('Number of filtered reads:'), '{:>10}'.format(filtered_count), file = logfile)
print('{:25}'.format('Number of excluded reads:'), '{:>10}'.format(excluded_count), '(excluded reads are listed below)', file = logfile)

#Calculate GC content
GC = (float(count_g) + float(count_c)) / float(nuc_sum)
#save nucleotide counts to logfile
print('\nTotal nucleotide counts\nA:', count_a, '\nT:', count_t, '\nG:', count_g, '\nC:', count_c, '\nGC content:', '{:.2f}'.format(GC), file = logfile)


####close files####
infile.close()
outfile.close()
logfile.close()



