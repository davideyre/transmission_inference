#!/usr/bin/env python
import sys, os, random, math
from time import clock, time
from optparse import OptionParser
from random import randint
from Bio import SeqIO

def IsInt(s):
	"""Does string represent an integer?"""
	try: 
		int(s)
		return True
	except ValueError:
		return False



class Node:
	"""Tree nodes"""
	samples=set()
	child1=None
	child2=None
	branch1=0.0
	branch2=0.0
	parent=None
	branchp=0.0
	subtree=""
	coalTime=0.0
	like=[] #likelihood for pruning algorithm
	disState=[] # List of discrete state values a the node.
	hist1=[] #mutational history to child 1
	hist2=[] #mutational history to child 2
	like1=[] #likelihood for state right before coalescent events, ignoring the coalescent event itself 
	like2=[]
	cacheInd=0.1



#samplList=[[2,2,3.5,'s1'],..], sampleList=[[numSamples,Deme,timeFromPres,samplingPointName],...]
#speciationList=[[1,2,3.5],...], speciationList=[[movingDeme,destinationDeme,speciationTime],..]
def SimuMultiSpCoal(sampleList,speciationList,popSizes,bottleneck,seedNum,demeCodes,commDeme):
	#Time in the input lists and output tree is expressed in terms of 2N generations, so to be consistent with msms. - over-riden to be in terms of N generations
	"""Simulate multi-species coalescent according to list of sampling events, speciation events, population sizes, and a bottleneck parameter."""
	nPop=len(popSizes)
	random.seed(seedNum)
	#sort samples into time order, most recent first
	sampleListSort=sorted(sampleList, key=lambda t: t[2])
	sampleList=sampleListSort
	
	#sort speciation events (i.e. transmission events), most recent first
	speciationListSort=sorted(speciationList, key=lambda t: t[2])
	speciationList=speciationListSort
	
	if len(speciationList)!=(len(popSizes)-1):
		print("Error: there must be numHosts-1 transmission (speciation) events!")
		exit()
	if len(sampleList)<2:
		print("Error: there must be at least 2 sampling events!")
		exit()
	
	sampleInd=0 #index to iterate through time-sorted samples list
	specInd=0 #index for speciation event list
	samples=[] #list of number of samples in each host
	demesOpen=[] #list of hosts that are currently infected
	popList=[] #list of nodes in each host
	
	for i in range(len(popSizes)):
		demesOpen.append(1) #set all hosts to be infected at present
		samples.append(0) #set all hosts to have zero samples
		popList.append([]) 
	
	
	#start coalescent process
	totSam=sum(samples)+len(sampleList)-sampleInd
	
	totalTime=0.0
	samCount=0 #overall number of samples obtained across all sampling points
	
	while totSam>1 or sampleInd!=len(sampleList):
		if sum(samples)>1:
			newSam=MultiSpCoalescentStep(samples,popSizes,totalTime)
			coalTime=newSam[2] #get the time of the next coalescent event
		else:
			coalTime=float("inf")
		
		if sampleInd<len(sampleList):
			sampleTime=sampleList[sampleInd][2] #get the time of the next sampling event
		else:
			sampleTime=float("inf")
		
		if specInd<len(speciationList):
			specTime=speciationList[specInd][2] #get the time of the next speciation event
		else:
			specTime=float("inf")
		
		#sample event
		if sampleTime<=coalTime and sampleTime<=specTime:
			#print "sample at time "+str(sampleTime)
			samples[sampleList[sampleInd][1]]+=sampleList[sampleInd][0] #add the number of samples obtained at this time
			totalTime=sampleList[sampleInd][2] #update the time to the new time
			totSam=sum(samples)
			
			if demesOpen[sampleList[sampleInd][1]]==0:
				#check sampling has occurred from infected host (open deme)
				print("Error: deme "+str(sampleList[sampleInd][1])+" already closed!")
				exit()
		
			#initialize leaf nodes
			for leaf in range(sampleList[sampleInd][0]): 
				#at sampleInd is current sampling point, sampleList[sampleInd][0] is the number of samples obtained
				#store sample as a leaf node in tree
				samCount+=1
				newnode=Node()
				newnode.samples=set([str(samCount)])
				newnode.disState=[sampleList[sampleInd][1]+1]
				#newnode.subtree=sampleList[sampleInd][3]+"_"+str(leaf+1)+"[&host=\""+demeCodes[sampleList[sampleInd][1]]+"\"]"#"[&host=\""+demeCodes[sampleList[sampleInd][1]]+"\",height="+str(totalTime)+"]"
				newnode.subtree=sampleList[sampleInd][3]+"[&host=\""+demeCodes[sampleList[sampleInd][1]]+"\"]"#"[&host=\""+demeCodes[sampleList[sampleInd][1]]+"\",height="+str(totalTime)+"]"
				newnode.coalTime=totalTime
				popList[sampleList[sampleInd][1]].append(newnode) #add new node to list of nodes within this host
				
			sampleInd+=1 #proceed to next sample
				
		#coalescent event
		elif coalTime<=specTime:
			#print "coalescent at time "+str(coalTime)
			node1=popList[newSam[0]].pop(newSam[1][1])
			node2=popList[newSam[0]].pop(newSam[1][0])
			newNode=Coalesce(node1,node2,newSam[2],demeCodes[newSam[0]]) # return new internal node
			popList[newSam[0]].append(newNode) #add new node to nodes in this host
			samples[newSam[0]]-=1 #update number of samples in this host
			totalTime=newSam[2] #update time
			totSam=sum(samples) #update number of samples
		
		#speciation event
		else:
			#print "speciation at time "+str(specTime)
			#specInd - index for most recent speciation event
			fromDeme=speciationList[specInd][0]
			toDeme=speciationList[specInd][1]
			bottleneckEvent = bottleneck
			if toDeme == commDeme:
				bottleneckEvent = bottleneck ##can edit here to enforce very high bottleneck for community imports
			if demesOpen[fromDeme]==0:
				print("Error: deme "+str(fromDeme)+" already closed!")
				exit()
			demesOpen[fromDeme]=0 #close deme 
			totalTime=speciationList[specInd][2] #update time to speciation event time
			specInd+=1 #increment index
			time=0.0
			oneDemePopSizes=[popSizes[fromDeme]] #population size in from deme
			#coalescent events in bottleneck - create a dummy amount of time and run coalescent events
			while time<bottleneckEvent:
				oneDemeSamples=[samples[fromDeme]]
				if oneDemeSamples[0]<2:
					break
				newSam=MultiSpCoalescentStep(oneDemeSamples,oneDemePopSizes,time)
				if time>bottleneckEvent:
					break
				else:
					time=newSam[2]
					node1=popList[fromDeme].pop(newSam[1][1])
					node2=popList[fromDeme].pop(newSam[1][0])
					newNode=Coalesce(node1,node2,totalTime,demeCodes[fromDeme])
					popList[fromDeme].append(newNode)
					samples[fromDeme]-=1
					totSam=sum(samples)
			#move lineages in closing deme to the destination deme
			samples[toDeme]+=samples[fromDeme]
			for i in range(samples[fromDeme]):
				popList[fromDeme][i].disState=[toDeme+1]
			for i in range(samples[fromDeme]):
				popList[toDeme].append(popList[fromDeme].pop())
			samples[fromDeme]=0

	#find the only node left and add ; to the tree
	pop=0
	while len(popList[pop])==0:
		pop+=1
	root=popList[pop][0]
	tree=root.subtree+";"
	return tree



#make a coalscent step for the multispecies coalescent
def MultiSpCoalescentStep(samples,popSizes,totalTime): 
	#coalescent rates
	coalRates=[]
	totalRate=0.0
	#get a coalescence rate for each host - n_samples * (n_samples - 1) 
	for pop in range(len(samples)):
		rate=samples[pop]*(samples[pop]-1)/(2*popSizes[pop])
		coalRates.append(rate)
		totalRate+=rate
	if totalRate<=0.0:
		time=float("inf")
		popToCoal=-1
		elemToCoal=[-1,-1]
	else:
		#sample time to next event
		time=random.expovariate(totalRate)
		#sample host to experience event
		popToCoal=SampleFromList(coalRates)
		#sample the 2 samples to coalesce
		elemToCoal=random.sample(range(samples[popToCoal]),2)
		elemToCoal.sort()
	#print 'Host: %s has potential coalescent event between %s and %s at time %s'%(popToCoal, elemToCoal[0], elemToCoal[1], time+totalTime)
	return [popToCoal,elemToCoal,time+totalTime] 
	#host to experiencing coalescent event, pair of samples coalescing, time of event




def Coalesce(node1,node2,time,deme): #Coalescing two nodes - store relevant data in node classes
	"""Coalesce two nodes into one"""
	newnode=Node()
	newnode.child1=node1
	newnode.child2=node2
	newnode.branch1=time-node1.coalTime
	newnode.branch2=time-node2.coalTime
	newnode.coalTime=time
	node1.parent=newnode
	node2.parent=newnode
	node1.branchp=time-node1.coalTime
	node2.branchp=time-node2.coalTime
	node1.subtree=node1.subtree+":"+str(time-node1.coalTime)
	node2.subtree=node2.subtree+":"+str(time-node2.coalTime)
	newnode.subtree="("+node1.subtree+","+node2.subtree+")"+"[&host=\""+deme+"\"]"#+"[&host=\""+deme+"\",height="+str(time)+"]"
	newnode.samples =(node1.samples | node2.samples)
	if node1.disState[0]==node2.disState[0]:
		newnode.disState=node1.disState[:1]
	else :
		print "error: coalescence of nodes with differing states."
	return newnode



def SampleFromList(listProbs): # sampling an element from a weighted list
	"""sample element according to a list of weights"""
	partial=0.0
	for element in range(len(listProbs)):
		partial+=listProbs[element]
	newList=[]
	for element in range(len(listProbs)):
		newList.append(listProbs[element]/partial)
	ranNum=random.random()
	partial=0.0
	for element in range(len(newList)):
		partial+=newList[element]
		if ranNum<partial:
			return element
	return len(listProbs)-1



def treeFromString(treeString,nodep): #Transform string into a bifurcating Tree structure without state at internal nodes
	"""Transform string into a bifurcating Tree structure without state at internal nodes (for now)."""
	treeString=treeString.replace(";",":0.0")
	if treeString[0]!="(":
		#Leaf case
		node=getLeaf(treeString)
		branch=getBranch(treeString)
		node.branchp=branch
		node.parent=nodep
		for sam in node.samples:
			break
		node.subtree=sam+":"+str(branch)
		return node	
	else:
		#internal node case
		indexEnd=len(treeString)-1
		while treeString[indexEnd]!=")":
			indexEnd-=1
			if indexEnd<0:
				print "parenthesis error while parsing tree"
				return None
		branch=getBranch(treeString[indexEnd+1:])
		node=Node()
		node.parent=nodep
		node.branchp=branch
		if treeString[indexEnd+1]=="[":
			#node.disState=[getState(treeString[indexEnd+1:])]
			ss=getState(treeString[indexEnd+1:])
			if IsInt(ss):
				node.disState=[int(ss)]
			else:
				node.disState=[ss]
		index=0
		countPar=0
		while treeString[index]!="," or countPar!=0:
			index+=1
			if treeString[index]=="(" or treeString[index]=="[":
				countPar+=1
			elif treeString[index]==")" or treeString[index]=="]":
				countPar-=1
		#Call iteratively function on subtrees	
		#print "subtree1", treeString[1:index]
		subtree1=treeFromString(treeString[1:index],node)
		node.child1=subtree1
		node.branch1=subtree1.branchp
		#print "subtree2", treeString[index+1:indexEnd]
		subtree2=treeFromString(treeString[index+1:indexEnd],node)
		node.child2=subtree2
		node.branch2=subtree2.branchp
		node.samples =(node.child1.samples | node.child2.samples)
		#print "New node",node.samples
		node.subtree= "(" + subtree1.subtree + "," + subtree2.subtree + "):" + str(branch)
		node.coalTime=node.child1.coalTime+node.branch1
		return node
		


def getLeaf(treeString): #Define node associated to tree leaf.
	"""Define node associated to tree leaf."""
	index=0
	while treeString[index]!=":":
		index+=1
	node=Node()
	if treeString[index-1]!="]":
		node.samples=set([treeString[:index]])
	else:
		index2=0
		while treeString[index2]!="[":
			index2+=1
		node.samples=set([treeString[:index2]])
		node.disState=[treeString[index2+1:index-1]]
		ss=treeString[index2+1:index-1]
		if IsInt(ss):
			node.disState=[int(ss)]
		else:
			node.disState=[ss]
	node.coalTime=0.0
	#print "New leaf", node.samples
	return node



def getState(treeString): #Get discrete state from string
	"""Get discret state from string"""
	index=0
	while treeString[index]!="]":
		index+=1
	index2=0
	while treeString[index2]!="[":
		index2+=1
	return treeString[index2+1:index]



def getBranch(treeString): #Get branch length from string
	"""Get branch length from string"""
	index=0
	while treeString[index]!=":":
		index+=1
	branch=float(treeString[index+1:])
	#print "branch length", branch
	return branch


def scaleTree(tree,scale,PNode):#scale branch lengths of a tree
	"""scale branch lengths of a tree."""
	if tree==None:
		return None
	newTree=Node()
	newTree.samples=tree.samples
	newTree.branch1=tree.branch1*scale
	newTree.branch2=tree.branch2*scale
	newTree.parent=PNode
	newTree.branchp=tree.branchp*scale
	newTree.coalTime=tree.coalTime*scale
	newTree.disState=tree.disState[:]
	hist=tree.hist1[:]
	for h in range(len(hist)):
		hist[h]=hist[h]*scale
	newTree.hist1=hist
	hist=tree.hist2[:]
	for h in range(len(hist)):
		hist[h]=hist[h]*scale
	newTree.hist2=hist
	
	newTree.child1=scaleTree(tree.child1,scale,newTree)
	newTree.child2=scaleTree(tree.child2,scale,newTree)
	if newTree.child1!=None:
		newTree.subtree="("+newTree.child1.subtree+":"+str(newTree.coalTime-newTree.child1.coalTime)+","+newTree.child2.subtree+":"+str(newTree.coalTime-newTree.child2.coalTime)+")"
		if newTree.disState!=[]:
			newTree.subtree+=("["+str(newTree.disState[0])+"]")
	else:
		for e in newTree.samples:
   			break
		newTree.subtree=str(e)+"["+str(newTree.disState[0])+"]"
	return newTree
	

def RemoveStates(tree):#Remove states from a tree string
	"""Remove the discrete states information (expressed as [etc.]) from a tree."""
	pos=0
	newTree=""
	while pos<len(tree):
		if tree[pos]=="[":
			while tree[pos]!="]":
				pos+=1
		else :
			newTree+=tree[pos]
		pos+=1
	return newTree

def getNexus(tree, sampleList):
	taxa = []
	for s in sampleList:
		for i in range(s[0]):
			taxa.append(s[3]+'_'+str(i+1))
	
	treeStringNexus = '#NEXUS\n'
	treeStringNexus += '\tDimensions ntax=%s;\n'%len(taxa)
	treeStringNexus += '\tTaxlabels\n';
	treeStringNexus += '\t\t%s'%'\n\t\t'.join(taxa)
	treeStringNexus += '\t\t;\n'
	treeStringNexus += 'End;\n\n'
	treeStringNexus += 'Begin trees;\n'
	treeStringNexus += 'tree TREE1 = %s\n'%tree
	treeStringNexus += 'End;\n'
	#add Nicola's formatting
	treeStringNexus +=  "begin figtree;\n	set appearance.branchColorAttribute=\"host\";\n	set appearance.branchColorGradient=true;\n	set nodeLabels.colorAttribute=\"host\";\n	set nodeLabels.displayAttribute=\"host\";\n	set nodeLabels.fontSize=13;\n 	set tipLabels.colorAttribute=\"host\";\n	set tipLabels.displayAttribute=\"host\";\n	set tipLabels.fontSize=13;\n end;\n"
	
	return treeStringNexus


if __name__ == "__main__":
	
	parser = OptionParser()
	parser.add_option( '-p', '--host_popsize', action = 'store', type='float', dest = 'hostPopSize', default = '1.0' )
	parser.add_option( '-c', '--community_popsize', action = 'store', type='float', dest = 'communityPopSize', default = '1000.0' )
	parser.add_option( '-b', '--bottleneck', action = 'store', type='float', dest = 'bottleneck', default = '1000.0' )
	parser.add_option( '-s', '--seed', action = 'store', type='int', dest = 'seed', default = 0 )
	parser.add_option( '-a', '--aln', action = 'store', type='int', dest = 'alnLength', default = 1000 )
	parser.add_option( '-m', '--mut', action = 'store', type='int', dest = 'mutationRate', default = 2 )
	parser.add_option( '-f', '--filepath', action = 'store', type='string', dest = 'pathDir', default = '.' )
	
	opts, args = parser.parse_args()
	if opts.seed == 0:
		opts.seed = randint(0,1e10)
		sys.stdout.write("Using seed %s\n"%opts.seed)
	
	#samplList=[[2,2,3.5,'s1'],..], sampleList=[[numSamples,Deme,timeFromPres,samplingPointName],...]
	#speciationList=[[1,2,3.5],...], speciationList=[[movingDeme,destinationDeme,speciationTime],..]
	#popSizes=[1.0,1.0,1.0,1.0]
	#sampleList=[[2,0,3.5, 's1'],[3,1,4.5, 's2'],[1,2,2.5, 's3'],[3,3,1.5, 's4']]
	#speciationList=[[0,1,4.6],[1,2,4.6],[2,3,24.6]]
	#bottleneck=1000.0
	#demeCodes = ['A', 'B', 'C', 'D']
	
	pathDir = opts.pathDir
	
	f = open('%s/siLog.txt'%pathDir, 'r')
	bin = f.next() #"patientid" "t_inf" "source" "source_type" "t_sample"
	sirLog = [l.strip().split() for l in f]
	f.close()
	
	#get maximum date - will be a sample time
	maxDate = max([int(i[4]) for i in sirLog])
	
	#encode patients as integers
	patientCode = dict()
	
	for i, patient in enumerate(sirLog):
		patientCode[patient[0]] = i
	
	nPatients = len(patientCode.keys())
	
	patientCode['-1'] = nPatients
	
	sampleList = [] #samplList=[[2,2,3.5,'s1'],..], sampleList=[[numSamples,Deme,timeFromPres,samplingPointName],...]
	demeCodes = []
	speciationList = [] #speciationList=[[1,2,3.5],...], speciationList=[[movingDeme,destinationDeme,speciationTime],..]
	popSizes = []
	
	for i in sirLog:
		sample = [1, patientCode[i[0]], maxDate-int(i[4]), "Samp_"+i[0]]
		sampleList.append(sample)
		
		demeCodes.append("Patient_"+i[0])
		
		speciation = [patientCode[i[0]], patientCode[i[2]], maxDate-int(i[1])]
		speciationList.append(speciation)
		popSizes.append(opts.hostPopSize)
	
	print sampleList
	demeCodes.append('Community')
	popSizes.append(opts.communityPopSize)
	#print popSizes
	#print demeCodes
	commDeme = len(demeCodes)-1

	tree=SimuMultiSpCoal(sampleList,speciationList,popSizes,opts.bottleneck,opts.seed,demeCodes,commDeme)
	treeString=tree
	treeStringNewick=RemoveStates(treeString)
	treeStringNexus = getNexus(treeString, sampleList)
	
	tree=treeFromString(tree,None)
	#tree=scaleTree(tree,2.0,None)
	print '\nHeight: %0.4f\n'%tree.coalTime #print height of tree
	
	#write trees in newick and nexus format
	f=open('%s/sim_newick.tree'%pathDir, 'w')
	f.write(treeStringNewick)
	f.close()
	
	f=open('%s/sim_nexus.tree'%pathDir, 'w')
	f.write(treeStringNexus)
	f.close()
	#print '%s\n'%treeString
	
	#run seq-gen - assume 1000 bases, at 2 base per year, rate 0.002 per year, per day 0.000005475701574
	mu = float(opts.mutationRate) / opts.alnLength / 365.25
	cmd = '/Applications/Seq-Gen.v1.3.3/seq-gen -mHKY -l%s -op -z%s -s%s %s/sim_newick.tree > %s/sim_alignment.phy'%(opts.alnLength, opts.seed, mu, pathDir, pathDir)
	sys.stdout.write(cmd+"\n")
	os.system(cmd+"\n")
	
	#convert output format to fasta
	aln = '%s/sim_alignment.phy'%pathDir
	aln_fa = '%s/sim_alignment.fa'%pathDir
	alnData = SeqIO.parse(aln, 'phylip')
	SeqIO.write(alnData, aln_fa, 'fasta')
		
	exit()

	