#!/usr/bin/env python3
import sys
import numpy as np

# Function parsing the PDB file
def parse_pdb(filename):
	pdb_coords_all={}
	pdb_coords_het={}
	f=open(filename,'r')
	for line in f:
		line=line.rstrip()
		# Check the field HETATM at the beginning of the line
		if line[:6]=='HETATM':
			chain=line[21]
			resn=int(line[22:27].strip())
			atom=line[13:17].strip()
			group=line[17:21].strip()
			x=float(line[30:38])
			y=float(line[38:46])
			z=float(line[46:54])
			coord=[x,y,z]
			#to have just the OXY and HEM groups 
			if group!='HOH':
			# Initialize the dictionary with heteroatoms coordinates of residue resn
				pdb_coords_het[resn]=pdb_coords_het.get(resn,{})
			# Add the atom's coordinates,the group and chain
				pdb_coords_het[resn][atom]=coord,group,chain
		elif line[:4]=='ATOM':
			chain=line[21]
			resn=int(line[22:26].strip())
			atom=line[12:16].strip()
			group=line[17:20].strip()
			x=float(line[30:38])
			y=float(line[38:46])
			z=float(line[46:54])
			coord=[x,y,z]
			# Initialize the dictionary with atoms coordinates of residue resn
			pdb_coords_all[resn]=pdb_coords_all.get(resn,{})
			# Add the atom's coordinates,the group and chain
			pdb_coords_all[resn][atom]=coord,group,chain
	#print(pdb_coords_het)
	return pdb_coords_all,pdb_coords_het


# Function calculating the distance between 2 points
def get_distance(coord1,coord2):
	return np.sqrt((coord1[0]-coord2[0])**2+\
                 (coord1[1]-coord2[1])**2+\
                 (coord1[2]-coord2[2])**2)


# Function returning the distances between all heteroatoms with all the other atoms
def get_all_dist(pdb_coords_all, pdb_coords_het):
	keys_all=list(pdb_coords_all.keys())
	keys_all.sort()
	keys_het=list(pdb_coords_het.keys())
	keys_het.sort()
	distances=[]
	for x in keys_all:
		for y in keys_het:
			for i in pdb_coords_all[x]:	#for each atom of the residue x
				for j in pdb_coords_het[y]:	#for each atom of the residue y
					dist=get_distance(pdb_coords_all[x][i][0],pdb_coords_het[y][j][0])
					if dist<=3.5 and dist!=0:
						dists=[]
						dists=[pdb_coords_all[x][i][2],x,pdb_coords_all[x][i][1],i,pdb_coords_het[y][j][2],y,pdb_coords_het[y][j][1],j,dist]
						distances.append(dists)
	new_dist=[]
	for i in distances:
		if i[3][0]=='C' or i[7][0]=='C': continue	#filter out the carbons
		
		if i[3][0]=='O' or i[3][0]=='N':	#to have electronegative atoms
			if i[7][0]=='O' or i[7][0]=='N':
				new_dist.append(i)
	df=pd.DataFrame(new_dist)
	df.to_csv('distances.tsv',sep='\t')
	#to get the unique residues
	uniqdist=[new_dist[0]]
	for i in range(1,len(new_dist)):
		if new_dist[i][1]==new_dist[i-1][1]:continue
		uniqdist.append(new_dist[i])
	dfuniq=pd.DataFrame(uniqdist)
	dfuniq.to_csv('unique_distances.tsv',sep='\t')
	return new_dist

def get_inter(pdb_coords_all):
	keys_all=list(pdb_coords_all.keys())
	keys_all.sort()
	distances=[]
	for x in keys_all:
		for y in keys_all:
			for i in pdb_coords_all[x]:	#for each atom of the residue x
				for j in pdb_coords_all[y]:	#for each atom of the residue y
					dist=get_distance(pdb_coords_all[x][i][0],pdb_coords_all[y][j][0])
					if dist<=4.5 and dist!=0:
						dists=[]
						dists=[pdb_coords_all[x][i][2],x,pdb_coords_all[x][i][1],i,pdb_coords_all[y][j][2],y,pdb_coords_all[y][j][1],j,dist]
						distances.append(dists)
	new_dist=[]
	for i in distances:
		if i[0]==i[4]: continue		#to have just interactions between different monomers
		if i[3][0]=='C' or i[7][0]=='C': continue	#filter out the carbons
		if i[3][0]=='O' or i[3][0]=='N':	#to have electronegative atoms
			if i[7][0]=='O' or i[7][0]=='N':
				new_dist.append(i)
	df=pd.DataFrame(new_dist)
	df.to_csv('interactions.tsv',sep='\t')
	return new_dist

# Function parsing the DSSP file
def parse_dssp(filename):
	dssp_list=[]
	# Initialize state variable c=0
	c=0
	f=open(filename)
	for line in f:
		line=line.rstrip()
		# Check for the beginning of data
		# Change the state variable c=1
		if line.find('  #  RESIDUE')>-1:
			c=1
			continue
		# Check if the state variable is 1
		# Select the correct chain
		if c==0: continue
		# Read all the variables resn, aa, acc, phi, psi
		aa=line[13]
		if aa.islower():
			aa='C'
		if aa=='!': continue
		resn=int(line[5:10])
		ss=line[16]
		if ss==' ':
			ss='C'
		acc=float(line[34:38])
		phi=float(line[103:109])
		psi=float(line[109:115])
		aa_dssp=[resn,aa,ss,acc,phi,psi]
		dssp_list.append(aa_dssp)
	return dssp_list
		

# Get sequence and secondary structure
def get_ss(dssp_list):		
	seq=''
	ss=''
	for aa in dssp_list:
		ss=ss+aa[2]
		seq=seq+aa[1]
	return seq,ss
	

# Count residues with a given secondary structure
def count_ss(ss,ss_type):
	c=0
	for i in ss:
		if i==ss_type: c=c+1
	return c 
	

# Return the list of all phi and psi angles
def get_ss_angle(dssp_list,ss_type):
	angles=[]
	for aa in dssp_list:
		if aa[2]==ss_type: angles.append((aa[-2],aa[-1]))
	return angles


# Return the list of residues accessibility
def get_aa_acc(dssp_list,aa_type):
	accs=[]
	for aa in dssp_list:
		if aa[1]==aa_type: accs.append(aa[3])
	return accs


def get_chain_acc(dssp_list):
	accs=[]
	for aa in dssp_list:
		accs.append(aa[3])
	return accs



if __name__ == '__main__':
	import pandas as pd
	'''filename=sys.argv[1]
	pdb_coords=parse_pdb(filename)
	
	#to calculate residues interacting with heteroatoms
	dists=get_all_dist(pdb_coords[0],pdb_coords[1])
	
	#to obtain all interacting atoms between chains
	final_dists=get_inter(pdb_coords[0])'''
	
	#print ('Dist for chain ', chain, ':',np.mean(ca_dists),np.std(ca_dists))
	Norm_Acc= {"A" :106.0,  "B" :160.0,         # D or N
   "C" :135.0,  "D" :163.0,  "E" :194.0,
   "F" :197.0,  "G" : 84.0,  "H" :184.0,
   "I" :169.0,  "K" :205.0,  "L" :164.0,
   "M" :188.0,  "N" :157.0,  "P" :136.0,
   "Q" :198.0,  "R" :248.0,  "S" :130.0,
   "T" :142.0,  "V" :142.0,  "W" :227.0,
   "X" :180.0,         # undetermined (deliberate)
   "Y" :222.0,  "Z" :196.0}         # E or Q
  
	ABCD=sys.argv[1]
	ABC= sys.argv[2]
	ABD= sys.argv[3]
	ACD= sys.argv[4]
	BCD= sys.argv[5]
	
	dssp_list_ABCD=parse_dssp(ABCD)
	dssp_list_ABC=parse_dssp(ABC)
	dssp_list_ABD=parse_dssp(ABD)
	dssp_list_ACD=parse_dssp(ACD)
	dssp_list_BCD=parse_dssp(BCD)
	
	
	
	#./parse_pdb_project.py ../data/1gzxABCD.dssp 1gzxABC.dssp 1gzxABD.dssp  1gzxACD.dssp 1gzxBCD.dssp
	print(get_chain_acc(dssp_list_ABCD))
	seq,ss=get_ss(dssp_list)
	print ('>SEQ\n%s' %seq)
	print ('>SS\n%s' %ss)
	print ('H count:',count_ss(ss,'H'))
	print ('E count:',count_ss(ss,'E'))
	h_angles=get_ss_angle(dssp_list,'H')
	e_angles=get_ss_angle(dssp_list,'E')
	print ('H angles:',np.mean([i[0] for i in h_angles]),np.mean([i[1] for i in h_angles]))
	print ('E angles:',np.mean([i[0] for i in e_angles]),np.mean([i[1] for i in e_angles]))
	print ('K acc:',np.mean(get_aa_acc(dssp_list,'K'))/Norm_Acc['K'])
	print ('V acc:',np.mean(get_aa_acc(dssp_list,'V'))/Norm_Acc['V'])


   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   

	
