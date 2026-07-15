import glob, os, sys

# Make one-hot matrix per each protein
multifastafile1 = sys.argv[1]
dir_name=multifastafile1.split('.')[0]
os.system('mkdir '+dir_name)
sequence = {}
counter = 0
seq = ''
header = ''
with open(multifastafile1, 'r') as f:
	for line in f:
		line = line.strip()
		if line.startswith('>'):
			if seq:
				sequence[header] = seq
			header = line[1:]
			counter += 1
			seq = ''
		else:
			seq += line

sequence[header] = seq
for name in sorted(sequence.keys()):
	SeqHead = name.split(' ')
	SeqHead2 = name.split('|')
	allername = SeqHead[0]
	allernameX = SeqHead2[0]
	allername = name
	allername = allername.rstrip(')')
	allernameX = allernameX.rstrip(')')
	SequenceName = sequence[name]
	mySeq = list(SequenceName)
	ll = len(SequenceName)
	
	with open(dir_name+"/"+f"{allernameX}.NN", 'w') as FH:
		FH.write(f"# File: {allername} \n")
		FH.write("#A R N D C Q E G H I L K M F P S T W Y V \n")
		for i in range(len(mySeq)):
			Resid = mySeq[i]
			FH.write(' '.join(['1' if Resid == aa else '0' for aa in 'ARNDCQEGHILKMFPSTWYV']) + "\n")

# zero-padding
flist=glob.glob(dir_name+'/*.NN')
cc=0
for f in flist:
	#print (f)
	inf=open(f,'r')
	with open(f, 'r') as inf_count:
		wct = sum(1 for _ in inf_count)
	#print (wct)
	if wct<1004:
		cc+=1
		outf=open(f+'.zeropad','w')

		ml=1000
		c=1
		while 1:
			line=inf.readline()
			if not line:break

			outf.write(line)
			outf.flush()
			if line[0]!='#':
				c+=1

		for i in range(1001-c):

			T=''
			for a in range(20):
				T+='0 '
			outf.write(T[:-1]+'\n')
			outf.flush()


# Concatenate to one one-hot encode matrix
flist=glob.glob(dir_name+'/*zeropad')

outf=open(dir_name+'_idx.txt','w')
outf2=open(dir_name+'.txt','w')
c=1
for f in flist:
	#print (f)
	num=os.path.basename(f).split('.')[0]
	inf=open(f,'r')
	outf.write(str(c)+'\t'+num+'\n')
	outf.flush()
	while 1:
		line=inf.readline()
		if not line:break
		if line[0]!='#':
			outf2.write(line)
			outf2.flush()
		c+=1
	c=c-2
