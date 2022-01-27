import sys
import argparse
import subprocess
import os
import tempfile
import glob

def collapse(genomic_range='', corrected_reads=''):
	parser = argparse.ArgumentParser(description='flair-collapse parse options', \
		usage='python flair.py collapse -g genome.fa -q <query.psl>|<query.bed> \
		-r <reads.fq>/<reads.fa> [options]')
	parser.add_argument('collapse')
	required = parser.add_argument_group('required named arguments')
	if not corrected_reads:
		required.add_argument('-q', '--query', type=str, default='', required=True, \
			action='store', dest='q', help='bed or psl file of aligned/corrected reads')
	required.add_argument('-g', '--genome', action='store', dest='g', \
		type=str, required=True, help='FastA of reference genome')
	required.add_argument('-r', '--reads', action='store', dest='r', nargs='+', \
		type=str, required=True, help='FastA/FastQ files of raw reads')
	parser.add_argument('-f', '--gtf', default='', action='store', dest='f', \
		help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
	parser.add_argument('-m', '--minimap2', type=str, default='minimap2', \
		action='store', dest='m', help='path to minimap2 if not in $PATH')
	parser.add_argument('-t', '--threads', type=int, \
		action='store', dest='t', default=4, help='minimap2 number of threads (4)')
	parser.add_argument('-p', '--promoters', action='store', dest='p', default='', \
		help='promoter regions bed file to identify full-length reads')
	parser.add_argument('--3prime_regions', action='store', dest='threeprime', default='', \
		help='TES regions bed file to identify full-length reads')
	parser.add_argument('-b', '--bedtools', action='store', dest='b', default='bedtools', \
		help='bedtools executable path, provide if TSS/TES regions specified and bedtools is not in $PATH')
	parser.add_argument('-sam', '--samtools', action='store', dest='sam', default='samtools', \
		help='samtools executable path if not in $PATH')
	parser.add_argument('-w', '--end_window', default='100', action='store', dest='w', \
		help='window size for comparing TSS/TES (100)')
	parser.add_argument('-s', '--support', default='3', action='store', dest='s', \
		help='''minimum number of supporting reads for an isoform;
		if s < 1, it will be treated as a percentage of expression of the gene (3)''')
	parser.add_argument('--stringent', default=False, action='store_true', dest='stringent', \
		help='''specify if all supporting reads need to be full-length \
		(80%% coverage and spanning 25 bp of the first and last exons)''')
	parser.add_argument('-n', '--no_redundant', default='none', action='store', dest='n', \
		help='''For each unique splice junction chain, report options include:
		none--best TSSs/TESs chosen for each unique set of splice junctions;
		longest--single TSS/TES chosen to maximize length;
		best_only--single most supported TSS/TES used in conjunction chosen (none)''')
	parser.add_argument('-i', '--isoformtss', default=False, action='store_true', dest='i', \
		help='when specified, TSS/TES for each isoform will be determined from supporting reads \
		for individual isoforms (default: not specified, determined at the gene level)')
	parser.add_argument('--no_gtf_end_adjustment', default=False, action='store_true', \
		dest='no_end_adjustment', \
		help='''when specified, TSS/TES from the gtf provided with -f will not be used to adjust 
		isoform TSSs/TESs each isoform will be determined from supporting reads''')
	parser.add_argument('--max_ends', default=2, action='store', dest='max_ends', \
		help='maximum number of TSS/TES picked per isoform (2)')
	parser.add_argument('--trust_ends', default=False, action='store_true', dest='trust_ends', \
		help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('--filter', default='default', action='store', dest='filter', \
		help='''Report options include:
		nosubset--any isoforms that are a proper set of another isoform are removed;
		default--subset isoforms are removed based on support;
		comprehensive--default set + all subset isoforms;
		ginormous--comprehensive set + single exon subset isoforms''')
	parser.add_argument('--fusion_dist', default=0, type=int, action='store', dest='fusion_dist', \
			help='''minimium distance between separate read alignments on the same chromosome to be
			considered a fusion, otherwise no reads will be assumed to be fusions''')
	parser.add_argument('--quality', type=int, action='store', dest='quality', default=1, \
		help='minimum MAPQ of read assignment to an isoform (1)')
	parser.add_argument('--keep_intermediate', default=False, action='store_true', \
		dest='keep_intermediate', \
		help='''specify if intermediate and temporary files are to be kept for debugging.
		Intermediate files include: promoter-supported reads file,
		read assignments to firstpass isoforms''')
	parser.add_argument('--generate_map', default=False, action='store_true', dest='generate_map', \
		help='''specify this argument to generate a txt file of read-isoform assignments 
		note: only works if the quantification method is not using salmon (default: not specified)''')
	parser.add_argument('--mm2_args', action='store', dest='mm2_args', \
		type=str, default='', help='''additional minimap2 arguments when aligning reads first-pass transcripts; 
		separate args by commas, e.g. --mm2_args=-I8g,--MD ''')
	parser.add_argument('--quiet', default=False, action='store_true', dest='quiet', \
			help='''Suppress progress statements from being printed''')
	parser.add_argument('--range', default='', action='store', dest='range', \
		help='''interval for which to collapse isoforms for, formatted chromosome:coord1-coord2 or 
		tab-delimited; if a range is specified, then the aligned reads bam must be specified with -r 
		and the query must be a sorted, bgzip-ed bed file''')
	parser.add_argument('--salmon', type=str, action='store', dest='salmon', \
		default='', help='Path to salmon executable, specify if salmon quantification is desired')
	parser.add_argument('--temp_dir', default='', action='store', dest='temp_dir', \
		help='directory for temporary files. use "./" to indicate current directory (default: python tempfile directory)')
	parser.add_argument('--split', default=False, action='store_true', dest='split', \
        help='''split sam file of collapsed reads due memory limitations (relevant only if salmon is NOT used''')
	parser.add_argument('-o', '--output', default='flair.collapse', \
		action='store', dest='o', help='output file name base for FLAIR isoforms (default: flair.collapse)')
	args, unknown = parser.parse_known_args()
	if unknown and not args.quiet:
		sys.stderr.write('Collapse unrecognized arguments: {}\n'.format(' '.join(unknown)))

	if corrected_reads:
		args.q = corrected_reads

	# housekeeping stuff
	tempfile_dir = tempfile.NamedTemporaryFile().name
	tempfile_name = 'tmpvcto_vk7'+'.'
	if args.temp_dir == '':
		args.temp_dir = tempfile_dir+'/'
	if not os.path.isdir(args.temp_dir):  # make temporary directory
		if subprocess.call(['mkdir', args.temp_dir]):
			sys.stderr.write('Could not make temporary directory {}\n'.format(args.temp_dir))
			return 1
	if args.temp_dir[-1] != '/':
		args.temp_dir += '/'

	if args.m[-8:] != 'minimap2':
		if args.m[-1] == '/':
			args.m += 'minimap2'
		else:
			args.m += '/minimap2'

	args.t, args.quality = str(args.t), str(args.quality)
	args.o += '.'
	min_reads = float(args.s) if float(args.s) >= 1 else 3

	if args.fusion_dist:
		args.trust_ends = True

	if not os.path.exists(args.q):
		sys.stderr.write('Query file path does not exist\n')
		return 1
	if os.stat(args.q).st_size == 0:
		sys.stderr.write('Query file is empty\n')
		return 1
	if float(args.s) < 1 and not args.f:
		sys.stderr.write('Provide gtf for gene grouping if -s is percentage of total gene expression\n')
		return 1
	# separate out the read sequences and corrected reads corresponding to the specified range
	ext = '.'+args.q[-3:]  # query file extension (bed or psl)
	precollapse = args.q  # query file unchanged
	args.r = args.r[0].split(',') if ',' in args.r[0] else args.r  # read sequences
	for r in args.r:
		if not os.path.exists(args.q):
			sys.stderr.write('Check that read file {} exists\n'.format(r))
			return 1

	count_files, align_files = [], []

	alignout = args.temp_dir + tempfile_name +'firstpass.'

	# count the number of supporting reads for each first-pass isoform
	if args.salmon:  # use salmon to count
		if subprocess.call([args.sam, 'view', '-F', '4', '-h', '-S', alignout+'sam'], \
			stdout=open(alignout+'mapped.sam', 'w')):
			return 1
		subprocess.call(['mv', alignout+'mapped.sam', alignout+'sam'])
		subprocess.call([args.salmon, 'quant', '--ont', '-t', args.o+'firstpass.fa', '-o',  alignout+'salmon', \
			'-p', args.t, '-l', 'A', '-a', alignout+'sam'], stderr=open(alignout+'salmon_stderr.txt', 'w'))
		count_file = alignout+'salmon/quant.sf'
		align_files += [alignout+'sam', alignout+'salmon/quant.sf']
	else:
		args.quality = '0' if args.trust_ends else args.quality
		if args.quality != '0':
			subprocess.call([args.sam, 'view', '-q', args.quality, '-h', '-S', alignout+'sam'], \
				stdout=open(alignout+'q.sam', 'w'), stderr=open(alignout+'q.samtools_stderr', 'w'))
			align_files += [alignout+'sam']
		else:
			subprocess.call(['mv', alignout+'sam', alignout+'q.sam'])
		count_cmd = [sys.executable, path+'bin/count_sam_transcripts.py', '-t', args.t, '--quality', args.quality]
		if args.stringent:
			count_cmd += ['--stringent', '-i', args.o+'firstpass'+ext]
		if args.trust_ends:
			count_cmd += ['--trust_ends']
		if args.generate_map:
			count_cmd += ['--generate_map', args.o+'isoform.read.map.txt']
		if args.fusion_dist:
			count_cmd += ['--fusion_dist', str(args.fusion_dist), ]
		if args.split:
			#subprocess.call(['split', '-C', '30GB', '-d', alignout+'q.sam', alignout+'q.sam.'])
			p1 = subprocess.Popen('grep -v ^@ ' + alignout+'q.sam', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			cmd = 'split -d -C 30G --filter="{ ' + args.sam + ' view -H ' + alignout + 'q.sam; cat; } > $FILE" - ' + alignout + 'r.sam.' 
			print(cmd) 
			p2 = subprocess.Popen(cmd, shell=True, stdin=p1.stdout)
			p1.stdout.close()
			out, err = p2.communicate()

			for file in os.listdir(args.temp_dir):
				filename = os.fsdecode(file)
				if filename.startswith(tempfile_name +'firstpass.' +'r.sam.'):
					suffix = filename[-2:]
					print(suffix)
					iocmd = ['-s', args.temp_dir + filename, '-o', alignout + 'r.counts.' + suffix]
					if subprocess.call(count_cmd + iocmd):
						sys.stderr.write('Failed at counting step for isoform read support\n')
						return 1
					count_files+=[alignout + 'r.counts.' + suffix]
					align_files+=[filename]
		else:
			iocmd = ['-s', alignout + 'q.sam', '-o', alignout + 'q.counts']
			if subprocess.call(count_cmd + iocmd):
				sys.stderr.write('Failed at counting step for isoform read support\n')
				return 1
			count_files = [alignout+'q.counts']
			align_files += [alignout+'q.sam']
	print(count_files)
	if subprocess.call([sys.executable, path+'bin/combine_counts.py'] + count_files + [args.o+'firstpass.r.counts']):
		sys.stderr.write('Failed at combining counts for transcripts\n')
		return 1

	if not args.quiet: sys.stderr.write('Filtering isoforms by read coverage\n')
	match_count_cmd = [sys.executable, path+'bin/match_counts.py', args.o+'firstpass.r.counts', \
		args.o+'firstpass'+ext, str(min_reads), args.o+'isoforms'+ext]
	if args.generate_map:
		match_count_cmd += ['--generate_map', args.o+'isoform.read.map.txt']
	subprocess.call(match_count_cmd)
	
""" 	if not args.range:  # also write .fa and .gtf files
		subprocess.call([sys.executable, path+'bin/psl_to_sequence.py', args.o+'isoforms'+ext, \
			args.g, args.o+'isoforms.fa'])
		if args.f:
			subprocess.call([sys.executable, path+'bin/psl_to_gtf.py', args.o+'isoforms'+ext], \
				stdout=open(args.o+'isoforms.gtf', 'w')) """

path = '/'.join(os.path.realpath(__file__).split("/")[:-1])+'/'
if len(sys.argv) < 2:
	sys.stderr.write('usage: python flair.py <mode> --help \n')
	sys.stderr.write('modes: align, correct, collapse, quantify, diffExp, diffSplice\n')
	sys.stderr.write('Multiple modules can be run when specified using numbers, e.g.:\n')
	sys.stderr.write('python flair.py 1234 ...\n')
	sys.exit(1)
else:
	mode = sys.argv[1].lower()

aligned_reads, corrected_reads, isoforms, isoform_sequences, counts_matrix = [0]*5

if mode == 'collapse' or ('3' in mode and '3.5' not in mode):
	if corrected_reads:
		status = collapse(corrected_reads=corrected_reads)
	else:
		status = collapse()
	if status == 1:
		sys.exit(1)
	else:
		isoforms, isoform_sequences = status