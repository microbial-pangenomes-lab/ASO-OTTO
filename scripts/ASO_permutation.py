import itertools as it
import sys
import argparse
import os
from io import StringIO as sio

def get_options():
	description = "Get gene cluster specific k-mers from a set of bacterial genomes"
	parser = argparse.ArgumentParser(description=description)

	parser.add_argument("-m", "--mismatches", type = int,
						default = 2,
						help = "number of mismatches (default: %(default)d)")

	parser.add_argument("-a", "--ASO", type = str,
						help = "ASO you want to get permutations of")

	parser.add_argument("-o", "--out", type = str,
						help = "output directory")

	return parser.parse_args()

# hacky complement DNA function. should be replaced with module from Bio!
def complement_dna(sequence):
	complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	complement_sequence = ''.join(complement_dict.get(base, base) for base in sequence)
	return complement_sequence

def main(mismatches, ASOfile, output):
	mima = mismatches
	out = output
	indexcount = 0
	template2coding = {}
	pupy = ["T","C","A","G"]
	buff=sio()

# reas ASOs and convert them from N-C oriented PNAs to 5'-3' DNA (coding strand) 
	with open(f"{ASOfile}","r") as AFI:
		for header in AFI:
			if header.startswith(">"):
				ASO = AFI.readline().lstrip().rstrip()
				revcomp = complement_dna(ASO[::-1])
				template2coding[f"{indexcount}"] = (f"{ASO}",f"{revcomp}")
			indexcount+=1


	for idx , temcod in template2coding.items():
		# split ASO into basic ASO that does not change
		# and pre-/suffix that will be permuted
		revcomp = temcod[1]
		pref = revcomp[:mima]
		suff = revcomp[-mima:]
		baseASO = revcomp[mima:-mima]
		preperm = []
		sufperm = []
		permset = set()

		mimadict={}
		for mismat in range(0,mima+1):
			mimadict[f"{mismat}"]=[]

		# build suffix permutated ASO variants
		sufperm.append([pref,])
		sufperm.append([baseASO,])
		for x in range(mima):
			sufperm.append(pupy)

		# build prefix permutated ASO variants
		for x in range(mima):
			preperm.append(pupy)
		preperm.append([baseASO,])
		preperm.append([suff,])

		if mima > 0:
			for permu in it.product(*preperm):
				permset.add("".join(permu))
			for permu in it.product(*sufperm):
				permset.add("".join(permu))
		else:
			permset.add(revcomp)
			
		# add ASO permutations to the mismatch dictionary
		# depending on the mismatching nucleotide positions
		for ASOperm in permset:
			if ASOperm==revcomp:
				mimadict["0"].append(ASOperm)
			elif ASOperm[1] != revcomp[1] or ASOperm[-2] != revcomp[-2]:
				mimadict["2"].append(ASOperm)
			else:
				mimadict["1"].append(ASOperm)

		# build dictionary of extensions for ASOs
		# that are shorter than 12 nucleotides
		# (required for query)
		if len(revcomp) != 12:
			dif = 12 - len(revcomp)
			difdict={}
			for x in range(1,dif+1):
				difdict[f"{x}"]=["".join(p) for p in [*it.product("".join(pupy),repeat=x)]]

			# extend ASOs that are shorter than 12 nucleotides
			for leng in difdict.keys():
				for mismat in mimadict.keys():
					numkey=int(leng)
					if numkey==dif:
						extpermset=set()
						for extperm in it.product(*[difdict[leng], mimadict[mismat]]):
							joiextper="".join(extperm)
							buff.write(f">{mismat}:{joiextper}:{idx}\n")
							buff.write(f"{joiextper}\n")
						for extperm in it.product(*[mimadict[mismat],difdict[leng]]):
							joiextper="".join(extperm)
							buff.write(f">{mismat}:{joiextper}:{idx}\n")
							buff.write(f"{joiextper}\n")
					else:
						for extperm in it.product(*[difdict[leng], mimadict[mismat], difdict[f"{dif-numkey}"]]):
							joiextper="".join(extperm)
							buff.write(f">{mismat}:{joiextper}:{idx}\n")
							buff.write(f"{joiextper}\n")
		else:
			for mismat in mimadict.keys():
				for perm in mimadict[mismat]:
					buff.write(f">{mismat}:{perm}:{idx}\n")
					buff.write(f"{perm}\n")

		# write permutations to file in fasta format
		with open(f"{out}/permutation.fasta", "w") as perfi:
			perfi.write(buff.getvalue())

	# write path to permutation.fasta to file. used as input
	# for the kmer-DB query
	with open(f"{out}/permutation_fasta.txt", "w") as interfi:
		interfi.write(f"{out}/permutation.fasta")

	# keep file with the original index as well as sequence 
	# in PNA-NC and DNA 5'-3' orientation
	with open(f"{out}/idx_template_coding.txt", "a+") as falist:
		falist.write(f"index,PNA_NC,DNA_coding\n")
		for idx , temcod in template2coding.items():
			falist.write(f"{idx},{temcod[0]},{temcod[1]}\n")

if __name__ == '__main__':
	args = get_options()
	out = args.out
	mima = args.mismatches
	ASOfile = args.ASO
	main(mima, ASOfile, out)