from Bio import AlignIO
import numpy as np
import os,sys,re,random,shutil,subprocess
from sklearn.linear_model import LinearRegression

def NULLELIM(cvs):
	aft_data = []
	for i in range(len(cvs)):
		temp = []
		temp = cvs[i]
		if re.search('.*\w+.*',temp):
			aft_data.append(temp)
		else:
			continue
	return aft_data

def FILEOPEN(data):
	seq = [];seq2 = [];seq_name = [];seq_name2 = []
	seq_data = (open(data).read()).split("\n\n")
	seq_data = NULLELIM(seq_data)
	num_of_data = len(seq_data)

	if num_of_data == int(1):
		seq_name = seq_data[0].split("\n")[0]
		seq_name = ''.join(re.split(">",seq_name)[1])
		seq = ''.join(seq_data[0].split("\n")[1:])
		seq_name2.append(seq_name)
		seq2.append(seq)
	else:
		for i in range(num_of_data):
			seq_name = seq_data[i].split("\n")[0]
			seq_name = ''.join(re.split(">",seq_name)[1])
			seq = ''.join(seq_data[i].split("\n")[1:])
			seq_name2.append(seq_name)
			seq2.append(seq)

	return seq_name2,seq2

def CCANALYSIS(standard_param, stp_seq, asr_seq, cc_list):
	stp_seq = list(stp_seq[0]); target_param = np.zeros(11,float)
	#print(stp_seq)
	stp_numA = stp_seq.count("A")
	stp_numC = stp_seq.count("C")
	stp_numE = stp_seq.count("E")
	stp_numH = stp_seq.count("H")
	stp_numI = stp_seq.count("I")
	stp_numK = stp_seq.count("K")
	stp_numQ = stp_seq.count("Q")
	stp_numT = stp_seq.count("T")
	stp_numV = stp_seq.count("V")
	stp_numW = stp_seq.count("W")
	stp_numY = stp_seq.count("Y")

	for i in range(len(asr_seq)):
		tmp_asrseq = []
		tmp_asrseq = list(asr_seq[i])
		tmp_asrnum = len(tmp_asrseq)

		asr_numA = tmp_asrseq.count("A"); asr_numC = tmp_asrseq.count("C"); asr_numE = tmp_asrseq.count("E")
		asr_numH = tmp_asrseq.count("H"); asr_numI = tmp_asrseq.count("I"); asr_numK = tmp_asrseq.count("K")
		asr_numQ = tmp_asrseq.count("Q"); asr_numT = tmp_asrseq.count("T"); asr_numV = tmp_asrseq.count("V")
		asr_numW = tmp_asrseq.count("W"); asr_numY = tmp_asrseq.count("Y")
		
		#print(asr_numA)
		#print(stp_numA)

		target_param[0] = ((asr_numA-stp_numA)/tmp_asrnum)*100.0
		target_param[1] = ((asr_numC-stp_numC)/tmp_asrnum)*100.0
		target_param[2] = ((asr_numE-stp_numE)/tmp_asrnum)*100.0
		target_param[3] = ((asr_numH-stp_numH)/tmp_asrnum)*100.0
		target_param[4] = ((asr_numI-stp_numI)/tmp_asrnum)*100.0
		target_param[5] = ((asr_numK-stp_numK)/tmp_asrnum)*100.0
		target_param[6] = ((asr_numQ-stp_numQ)/tmp_asrnum)*100.0
		target_param[7] = ((asr_numT-stp_numT)/tmp_asrnum)*100.0
		target_param[8] = ((asr_numV-stp_numV)/tmp_asrnum)*100.0
		target_param[9] = ((asr_numW-stp_numW)/tmp_asrnum)*100.0
		target_param[10] = ((asr_numY-stp_numY)/tmp_asrnum)*100.0

		#target_param = (target_param - np.average(target_param))/np.std(target_param)
		print(target_param)

		cc_score = np.corrcoef(target_param,standard_param)[0][1]
		cc_list[i] = cc_score

	return cc_list

while len(sys.argv)>1:
	option = sys.argv[1]
	del sys.argv[1]
	if option == '-STP':
		stp_data = sys.argv[1]
		del sys.argv[1]
	elif option == '-ASRSEQ':
		asr_data = sys.argv[1]
		del sys.argv[1]
	elif option == '-OUTPUT':
		output = sys.argv[1]
		del sys.argv[1]

stp_name, stp_seq = FILEOPEN(stp_data)
asr_name, asr_seq = FILEOPEN(asr_data)
cc_list = np.zeros(len(asr_seq))

#stp_name and stp_seq are formed by only one element, whereas the asr_name and asr_seq contain
#more than one of elements.
print(stp_name,stp_seq)
print(asr_name,asr_seq)

#Mean difference for num. of AAs per thousands in mesophiles and thermophiles (thermophiles-mesophiles)
#standard_param = [A, C, E, H, I, K, Q, T, V, W, Y] of which p-value were <0.05.
#Gregory A.C. Singer and Donal A. Hickey, Gene 317, 2003, 39-47
standard_param = np.array([-21.6, -3.1, 22.6, -7.5, 16.3, 20.9, -25.3, -9.5, 11.1, -2.7, 9.4])
#standard_param = (standard_param-np.average(standard_param))/np.std(standard_param)
print(standard_param)

cc_list = CCANALYSIS(standard_param, stp_seq, asr_seq, cc_list)

print(asr_name)
print(cc_list)

#Data output
outputfile = open(output,'w')
outputfile.write("#The correlation of the prameters between standard and target protein.\n")
outputfile.write("#Mean difference for num. of AAs per thousands in mesophiles and thermophiles (thermophiles-mesophiles)\n")
outputfile.write("#standard_param = [A, C, E, H, I, K, Q, T, V, W, Y] of which p-value were <0.05.\n")
outputfile.write("#ref. Gregory A.C. Singer and Donal A. Hickey, Gene 317, 2003, 39-47\n\n")
outputfile.write("#STP name is as follows:\n")
outputfile.write(f">{stp_name[0]}\n{stp_seq[0]}\n\n")

num_of_list = len(cc_list)

for i in range(num_of_list):
	tmp_argmax = np.argmax(cc_list)
	tmp_cc = cc_list[tmp_argmax]
	tmp_seq = asr_seq[tmp_argmax]
	tmp_seqname = asr_name[tmp_argmax]
	cc_list[tmp_argmax] = -100000

	outputfile.write(f"#Rank {i+1}: CC_value was {tmp_cc}\n")
	outputfile.write(f">{tmp_seqname}\n{tmp_seq}\n\n")












