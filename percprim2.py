import sys
import random

class errDone(Exception):
	pass

# http://physics.nist.gov/cuu/Constants/index.html
m_n = 1.0086649
m_p = 1.007276
m_e = 0.00054857991

# -------------------- BINDING ENERGY PER NUCLEON ------------------------ #

# http://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl
# http://amdc.in2p3.fr/masstables/Ame2003/mass.mas03round
# http://amdc.in2p3.fr/web/masseval.html
fo=open("masses.txt","r")
lines=fo.readlines()
fo.close()
lines=lines[1:]	# Get rid of header

Amass=[100000000]*120
Aappr=[1]*120
for i in range(120):
	Amass[i]=[100000000]*300
	Aappr[i]=[1]*300
# Amass[Z][N] = mass(N,Z)
# Aappr[Z][N] = 1 for approximation; 0 for experimentally accurate

for line in lines:
	vals=line.split("\t")	# 
	vals[-1]=vals[-1][:-1]	# Get rid of '\n'
	N=int(vals[0])
	Z=int(vals[1])
	if vals[3][-1]=="#":
		Amass[Z][N]=(float(vals[2]+"."+vals[3][:-1]))
		Aappr[Z][N]=1
	else:
		Amass[Z][N]=(float(vals[2]+"."+vals[3]))
		Aappr[Z][N]=0

def B(N,Z):
	if N==0 and Z==1:
		return 0
	if N==1 and Z==0:
		return +0.00084
	if N>299 or Z>119:
		return 100000000
	if N<0 or Z<0:
		return 100000000
	if N==0 and Z==0:
		return 100000000
	mf=Amass[Z][N]
	mi=m_n*N + m_p*Z + m_e*Z
	return float(mf-mi)/(N+Z)

print "Calculating... Please wait."
print "Press CTRL+C at any time to pause running and view current composition of the universe"
print "If you do so, you will be able to resume calculation by typing 'y' or 'Y'"
#print Amass
# --------------------------- NAME OF ELEMENT --------------------------- #

fo=open("ele.txt")
lines=fo.readlines()
fo.close()

Ele={0:"n"}
for line in lines:
	vals=line.split("\t")
	Ele[int(vals[0])]=vals[1]

# ------------------------------------------------------------------------ #

def do():
	print
	print min([min(k) for k in Nctr])
	tot=0
	eles=[]
	eless=[]
	for j in range(120):
		for i in range(300):
			tot+=Nctr[i][j]
			if Nctr[i][j]>0:
				eles.append([Nctr[i][j],i,j])
			if Nmade[i][j]==1:
				eless.append([Nctr[i][j],i,j])
	
	print "Elements made at some time or another:"
	for elem in eless:
		print str(elem[1]+elem[2])+Ele[elem[2]],
	
	print
	print
	print "Current percentage composition:"
	eles.sort()
	eles.reverse()
	for elem in eles:
		print "%.2f%% %d%s; " % (float(elem[0]*100.0/tot),elem[1]+elem[2],Ele[elem[2]]) ,
	print

# ------------------------------------------------------------------------ #

Nctr=[0]*300
Nmade=[0]*300
for i in range(300):
	Nctr[i]=[0]*120
	Nmade[i]=[0]*120
# Nctr[n][z]

Ne=[[1,0]]*5000+[[0,1]]*5000	# [N,Z]
Nctr[1][0]+=5000
Nctr[0][1]+=5000
Nmade[1][0]=1
Nmade[0][1]=1
#print Ne

random.seed()
ctr=0
while(1):
	try:
		pos=random.randint(0,len(Ne)-1)
		N1=Ne[pos][0]
		Z1=Ne[pos][1]
		
		pos2=random.randint(0,len(Ne)-1)
		N2=Ne[pos2][0]
		Z2=Ne[pos2][1]
		
		if pos==pos2:
			continue
		
		plausible=[]	# Using the selected atoms, what all can be made
		Bplaus=[]		# Their binding energies
		Wplaus=[]
		
		# Beta-
		if B(N1-1,Z1+1) <= B(N1,Z1):
			plausible.append([N1-1,Z1+1])
			Bplaus.append(B(N1-1,Z1+1))
			Wplaus.append(1)
		# Beta+
		if B(N1+1,Z1-1) <= B(N1,Z1):
			plausible.append([N1+1,Z1-1])
			Bplaus.append(B(N1+1,Z1-1))
			Wplaus.append(1)
		# Alpha
		if B(N1-2,Z1-2) <= B(N1,Z1):
			plausible.append([N1-2,Z1-2])
			Bplaus.append(B(N1-2,Z1-2))
			Wplaus.append(1)
		
		# Beta-
		if B(N2-1,Z2+1) <= B(N2,Z2):
			plausible.append([N2-1,Z2+1])
			Bplaus.append(B(N2-1,Z2+1))
			Wplaus.append(2)
		# Beta+
		if B(N2+1,Z2-1) <= B(N2,Z2):
			plausible.append([N2+1,Z2-1])
			Bplaus.append(B(N2+1,Z2-1))
			Wplaus.append(2)
		# Alpha
		if B(N2-2,Z2-2) <= B(N2,Z2):
			plausible.append([N2-2,Z2-2])
			Bplaus.append(B(N2-2,Z2-2))
			Wplaus.append(2)
		
		# Fusion - No biproducts
		if B(N1+N2,Z1+Z2) <= B(N1,Z1) and B(N1+N2,Z1+Z2) <= B(N2,Z2):
			plausible.append([N1+N2,Z1+Z2])
			Bplaus.append(B(N1+N2,Z1+Z2))
			Wplaus.append(0)
		
		ind=-1
		if len(plausible)>0:
			ind=Bplaus.index(min(Bplaus))
			#print pos, pos2,Bplaus[ind]
			if Bplaus[ind]>=0:
				pass
			elif Wplaus[ind]==0:
				Ne.pop(Ne.index([N1,Z1]))
				Nctr[N1][Z1]-=1
				Ne.pop(Ne.index([N2,Z2]))
				Nctr[N2][Z2]-=1
				Ne.append(plausible[ind])
				Nctr[plausible[ind][0]][plausible[ind][1]]+=1
				Nmade[plausible[ind][0]][plausible[ind][1]]=1
				#print "Z=",Z1,"; N=",N1,"; A=",Z1+N1,";R",Ele[Z1]
				#print "Z=",Z2,"; N=",N2,"; A=",Z2+N2,";R",Ele[Z2]
				#print "Z=",plausible[ind][1],"; N=",plausible[ind][0],"; A=",plausible[ind][1]+plausible[ind][0],";P",Ele[plausible[ind][1]]
				#continue
			elif Wplaus[ind]==1:
				Ne[pos]=plausible[ind]
				Nctr[N1][Z1]-=1
				Nctr[plausible[ind][0]][plausible[ind][1]]+=1
				Nmade[plausible[ind][0]][plausible[ind][1]]=1
				#print "Z=",Z1,"; N=",N1,"; A=",Z1+N1,";R",Ele[Z1]
				#print "Z=",plausible[ind][1],"; N=",plausible[ind][0],"; A=",plausible[ind][1]+plausible[ind][0],";P",Ele[plausible[ind][1]]
				#continue
			elif Wplaus[ind]==2:
				Ne[pos2]=plausible[ind]
				Nctr[N2][Z2]-=1
				Nctr[plausible[ind][0]][plausible[ind][1]]+=1
				Nmade[plausible[ind][0]][plausible[ind][1]]=1
				#print "Z=",Z2,"; N=",N2,"; A=",Z2+N2,";R",Ele[Z2]
				#print "Z=",plausible[ind][1],"; N=",plausible[ind][0],"; A=",plausible[ind][1]+plausible[ind][0],";P",Ele[plausible[ind][1]]
				#continue
			else:
				#print N1,Z1,N2,Z2
				raise Exception
		else:
			pass
			#print "Z=",Z1,"; N=",N1,"; A=",Z1+N1,";R",Ele[Z1]
			#print "Z=",Z2,"; N=",N2,"; A=",Z2+N2,";R",Ele[Z2]
		#print "--"
		if ind!=-1:
			if Bplaus[ind]==100000000:
				break
		
		if ind==-1:
			ctr+=1
		else:
			ctr=0
		
		if ctr==5000:
			raise errDone
		'''
		if ind!=-1:
			if plausible[ind]==[3,2]:
				raw_input()
		'''
	except KeyboardInterrupt:
		print
		print ctr, "turns since last reaction."
		do()
		a=raw_input("Go on?")
		if len(a)==0:
			break
		elif a[0]=='y' or a[0]=='Y':
			pass
		else:
			break
	except errDone:
		print "Counter reached maximum value."
		print "ctr =",ctr
		print
		break

do()
#print Ne