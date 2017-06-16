import numpy as np
import csv
import matplotlib.pyplot as plt


popsize = ['10','20','30','40','50']
safesites = ['0.25','0.3','0.35','0.4','0.45','0.5']
#safesites = ['0.25','0.3','0.35']
probnews = ['0.0','0.2','0.4','0.6','0.8','1.0']
#probnews = ['0.0','0.2','0.4','0.6']


# demographic files 

persist = np.zeros((len(popsize),len(safesites),len(probnews)))
counti = 0
matrices = []        
for i in popsize:
	matj = []
	countj = 0
	for j in safesites:
		mati = []
		countk = 0
		for k in probnews:
			mat1 = []
			file1 = open('demo_' + i + j + k + '.dat', 'rb')
			reader1 = csv.reader(file1, delimiter='\t')
			for row in reader1:
				mat1.append(row)
			persistence = mat1[0:11]
			if persistence:
				persist[counti,countj,countk] = float(persistence[7][0].strip("Average Persistence: "))	
			del mat1[0:11]
			mat11 =[]
			for l in mat1:
				if len(l) > 1:
					mat11.append(l)        
			npmat1 = np.zeros((len(mat11),14))
			for ii in range(len(mat11)):
				for jj in range(14):
					if mat11[ii][jj] == '.':
						mat11[ii][jj] = 0
					npmat1[ii,jj] = float(mat11[ii][jj])
			file1.close() 
			mati.append(npmat1)
			countk += 1
		countj += 1
		matj.append(mati)
	matrices.append(matj)
	counti += 1

# fitness files

counti = 0
matricesfit = []        
for i in popsize:
	matj = []
	countj = 0
	for j in safesites:
		mati = []
		countk = 0
		for k in probnews:
			mat1 = []
			file1 = open('fit__' + i + j + k + '.dat', 'rb')
			reader1 = csv.reader(file1, delimiter='\t')
			for row in reader1:
				mat1.append(row)
			del mat1[0:11]
			mat11 =[]
			for l in mat1:
				if len(l) > 1:
					mat11.append(l)        
			npmat1 = np.zeros((len(mat11),17))
			for ii in range(len(mat11)):
				for jj in range(17):
					if mat11[ii][jj] == '.':
						mat11[ii][jj] = 0
					npmat1[ii,jj] = float(mat11[ii][jj])
			file1.close() 
			mati.append(npmat1)
			countk += 1
		countj += 1
		matj.append(mati)
	matricesfit.append(matj)
	counti += 1


# genetics files

counti = 0
matricesgen = []        
for i in popsize:
	matj = []
	countj = 0
	for j in safesites:
		mati = []
		countk = 0
		for k in probnews:
			mat1 = []
			file1 = open('gen__' + i + j + k + '.dat', 'rb')
			reader1 = csv.reader(file1, delimiter='\t')
			for row in reader1:
				mat1.append(row)
			del mat1[0:8]
			mat11 =[]
			for l in mat1:
				if len(l) > 1:
					mat11.append(l)        
			npmat1 = np.zeros((len(mat11),16))
			for ii in range(len(mat11)):
				for jj in range(16):
					if mat11[ii][jj] == '.':
						mat11[ii][jj] = 0
					npmat1[ii,jj] = float(mat11[ii][jj])
			file1.close() 
			mati.append(npmat1)
			countk += 1
		countj += 1
		matj.append(mati)
	matricesgen.append(matj)
	counti += 1

### Control simulations
matc = []
filec = open('demo_control2.dat', 'rb')
readerc = csv.reader(filec, delimiter='\t')
for row in readerc:
	matc.append(row)
 
del matc[0:11]
matcc =[]
for l in matc:
	if len(l) > 1:
		matcc.append(l)        
npmat_control = np.zeros((len(matcc),14))
for ii in range(len(matcc)):
	for jj in range(14):
		if matcc[ii][jj] == '.':
			matcc[ii][jj] = 0
		npmat_control[ii,jj] = float(matcc[ii][jj])
filec.close() 	

matc = []
filec = open('gen_control2.dat', 'rb')
readerc = csv.reader(filec, delimiter='\t')
for row in readerc:
	matc.append(row)
 
del matc[0:8]
matcc =[]
for l in matc:
	if len(l) > 1:
		matcc.append(l)        
npmat_control_g = np.zeros((len(matcc),16))
for ii in range(len(matcc)):
	for jj in range(16):
		if matcc[ii][jj] == '.':
			matcc[ii][jj] = 0
		npmat_control_g[ii,jj] = float(matcc[ii][jj])
filec.close() 

matc = []
filec = open('fit_control2.dat', 'rb')
readerc = csv.reader(filec, delimiter='\t')
for row in readerc:
	matc.append(row)
 
del matc[0:11]
matcc =[]
for l in matc:
	if len(l) > 1:
		matcc.append(l)        
npmat_control_f = np.zeros((len(matcc),16))
for ii in range(len(matcc)):
	for jj in range(16):
		if matcc[ii][jj] == '.':
			matcc[ii][jj] = 0
		npmat_control_f[ii,jj] = float(matcc[ii][jj])
filec.close() 


### --- Create a matrices of last value of temporal dynamics

# Demographics

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for k in range(len(probnews)):
	mateavail = np.zeros((len(popsize),len(safesites)))
	for i in range(len(popsize)):
		for j in range(len(safesites)):
			mateavail[i,j] = np.mean(matrices[i][j][k][-10:-1,12])
	axs[k].imshow(mateavail,interpolation="nearest", origin='lower');
	axs[k].set_title(str(probnews[k]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[k].imshow(mateavail,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()
		

fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	mateavail = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			mateavail[j,k] = np.mean(matrices[i][j][k][-1,12])
	axs[i].imshow(mateavail,interpolation="nearest", origin='lower');
	axs[i].set_title("Pop size = "+str(popsize[i]),fontsize=16)
	axs[i].set_xticks([0.0, 1.0,2.0,3.0,4.0,5.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	if i == 2:
		axs[i].set_xlabel('Prob new S alleles (GR)',fontsize=16)
	axs[i].set_ylabel('Safe sites increase (SR)',fontsize=16)
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
fig.suptitle("# Compatible reproductive individuals",fontsize=20)	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(mateavail,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()	

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for k in range(len(probnews)):
	reproinds = np.zeros((len(popsize),len(safesites)))
	for i in range(len(popsize)):
		for j in range(len(safesites)):
			reproinds[i,j] = np.mean(matrices[i][j][k][-10:-1,4])
	axs[k].imshow(reproinds,interpolation="nearest", origin='lower');
	axs[k].set_title(str(probnews[k]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[k].imshow(reproinds,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()


fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	reproinds = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			reproinds[j,k] = np.mean(matrices[i][j][k][-1,4])
	axs[i].imshow(reproinds,interpolation="nearest", origin='lower');
	axs[i].set_title("Pop size = "+str(popsize[i]),fontsize=16)
	axs[i].set_xticks([0.0, 1.0,2.0,3.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	if i == 2:
		axs[i].set_xlabel('Prob new S alleles (GR)',fontsize=16)
	axs[i].set_ylabel('Safe sites increase (SR)',fontsize=16)
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
		
fig.suptitle("Reproductive individuals",fontsize=20)	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(reproinds,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()	


fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	axs[i].imshow(persist[i],interpolation="nearest", origin='lower',cmap=plt.get_cmap('YlGn'));
	axs[i].set_title("Pop size = "+str(popsize[i]),fontsize=20)
	axs[i].set_xticks([0.0, 1.0,2.0,3.0,4.0,5.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	if i == 2:
		axs[i].set_xlabel('Prob new S alleles (GR)',fontsize=16)
	axs[i].set_ylabel('Safe sites increase (SR)',fontsize=16)
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
		
fig.suptitle("Persistence",fontsize=20)	
fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.025])
fig.colorbar(axs[0].imshow(persist[0],interpolation="nearest", origin='lower',cmap=plt.get_cmap('YlGn')), cax=cbar_ax,orientation='horizontal')
fig.tight_layout() 
plt.show()	

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for k in range(len(probnews)):
	reproinds = np.zeros((len(popsize),len(safesites)))
	for i in range(len(popsize)):
		for j in range(len(safesites)):
			reproinds[i,j] = np.mean(matrices[i][j][k][-10:-1,10])
	axs[k].imshow(reproinds,interpolation="nearest", origin='lower');
	axs[k].set_title(str(probnews[k]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[k].imshow(reproinds,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()
		
reproinds = np.zeros((len(popsize),len(safesites)))
for i in range(len(popsize)):
	for j in range(len(safesites)):
		reproinds[i,j] = np.mean(matrices[i][j][1][-10:-1,4])
plt.imshow(reproinds,interpolation="nearest", origin='lower'); plt.ylabel("Pop increase");plt.xlabel("Safe sites");plt.colorbar(); plt.show()

veginds = np.zeros((len(popsize),len(safesites)))
for i in range(len(popsize)):
	for j in range(len(safesites)):
		veginds[i,j] = np.mean(matrices[i][j][1][-10:-1,2])
plt.imshow(veginds,interpolation="nearest", origin='lower'); plt.ylabel("Pop increase");plt.xlabel("Safe sites");plt.colorbar(); plt.show()

lrep = np.zeros((len(popsize),len(safesites)))
for i in range(len(popsize)):
	for j in range(len(safesites)):
		lrep[i,j] = np.mean(matrices[i][j][1][-10:-1,10])
plt.imshow(lrep,interpolation="nearest", origin='lower'); plt.ylabel("Pop increase");plt.xlabel("Safe sites");plt.colorbar(); plt.show()
				
# Fitness				
sds1 = np.zeros((len(popsize),len(safesites)))
for i in range(len(popsize)):
	for j in range(len(safesites)):
		sds1[i,j] = np.mean(matricesfit[i][j][3][-10:-1,1])
plt.imshow(sds1,interpolation="nearest", origin='lower'); plt.ylabel("Pop increase");plt.xlabel("Safe sites");plt.colorbar(); plt.show()

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for k in range(len(probnews)):
	vsds1 = np.zeros((len(popsize),len(safesites)))
	for i in range(len(popsize)):
		for j in range(len(safesites)):
			vsds1[i,j] = np.mean(matricesfit[i][j][k][-10:-1,11])
	axs[k].imshow(vsds1,interpolation="nearest", origin='lower');
	axs[k].set_title(str(probnews[k]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[k].imshow(vsds1,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()	

fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	sds1 = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			sds1[j,k] = np.mean(matricesfit[i][j][k][-1,5])
	axs[i].imshow(sds1,interpolation="nearest", origin='lower');
	axs[i].set_title("Pop size = "+str(popsize[i]))
	axs[i].set_xticks([0.0, 1.0,2.0,3.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	axs[i].set_xlabel('Prob new S alleles (GR)')
	axs[i].set_ylabel('Safe sites increase (SR)')
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
		
fig.suptitle("Average seeds/plant before dispersal",fontsize=20)	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(sds1,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()

fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	sds1 = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			sds1[j,k] = np.mean(matricesfit[i][j][k][-1,11])
	axs[i].imshow(sds1,interpolation="nearest", origin='lower');
	axs[i].set_title("Pop size = "+str(popsize[i]),fontsize=16)
	axs[i].set_xticks([0.0, 1.0,2.0,3.0,4.0,5.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	if i == 2:
		axs[i].set_xlabel('Prob new S alleles (GR)',fontsize=16)
	axs[i].set_ylabel('Safe sites increase (SR)',fontsize=16)
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
		
fig.suptitle("Average seeds/plant post dispersal",fontsize=20)	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(sds1,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()		

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for k in range(len(probnews)):
	vsds1 = np.zeros((len(popsize),len(safesites)))
	for i in range(len(popsize)):
		for j in range(len(safesites)):
			vsds1[i,j] = np.mean(matricesfit[i][j][k][-10:-1,1])
	axs[k].imshow(vsds1,interpolation="nearest", origin='lower');
	axs[k].set_title(str(probnews[k]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[k].imshow(vsds1,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()	

# Genetics
sgene = np.zeros((len(popsize),len(safesites)))
for i in range(len(popsize)):
	for j in range(len(safesites)):
		sgene[i,j] = np.mean(matricesgen[i][j][1][-10:-1,8])
plt.imshow(sgene,interpolation="nearest", origin='lower'); plt.ylabel("Pop increase");plt.xlabel("Safe sites");plt.colorbar(); plt.show()

fig, axs = plt.subplots(2,3, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for k in range(len(probnews)):
	sgene = np.zeros((len(popsize),len(safesites)))
	for i in range(len(popsize)):
		for j in range(len(safesites)):
			sgene[i,j] = np.mean(matricesgen[i][j][k][-1,4])
	axs[k].imshow(sgene,interpolation="nearest", origin='lower');
	axs[k].set_title(str(probnews[k]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[k].imshow(sgene,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()			

fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	sgene = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			sgene[j,k] = np.mean(matricesgen[i][j][k][-1,14])
	axs[i].imshow(sgene,interpolation="nearest", origin='lower');
	axs[i].set_title("Pop size = "+str(popsize[i]))
	axs[i].set_xticks([0.0, 1.0,2.0,3.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	axs[i].set_xlabel('Prob new S alleles (GR)')
	axs[i].set_ylabel('Safe sites increase (SR)')
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
		
fig.suptitle("Fis",fontsize=20)	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(sgene,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()


fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	sgene = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			sgene[j,k] = np.mean(matricesgen[i][j][k][-1,10])
	axs[i].imshow(sgene,interpolation="nearest", origin='lower');
	axs[i].set_title("Pop size = "+str(popsize[i]))
	axs[i].set_xticks([0.0, 1.0,2.0,3.0,4.0,5.0])
	axs[i].set_xticklabels(probnews)
	axs[i].set_yticklabels(['0.0', '0.25','0.3','0.35','0.4','0.45','0.5'])
	if i == 2:
		axs[i].set_xlabel('Prob new S alleles (GR)',fontsize=14)
	axs[i].set_ylabel('Safe sites increase (SR)',fontsize=14)
	if i > 0:
		axs[i].get_yaxis().set_visible(False)
		
fig.suptitle("Average # S alleles",fontsize=20)	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(sgene,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()

fig, axs = plt.subplots(1,5, figsize=(15, 6), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .5, wspace=.001)

axs = axs.ravel()
for i in range(len(popsize)):
	sgene = np.zeros((len(safesites),len(probnews)))
	for j in range(len(safesites)):
		for k in range(len(probnews)):
			sgene[j,k] = np.mean(matricesgen[i][j][k][-1,2])
	axs[i].imshow(sgene,interpolation="nearest", origin='lower');
	axs[i].set_title(str(popsize[i]))
	
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(axs[i].imshow(sgene,interpolation="nearest", origin='lower'), cax=cbar_ax)
plt.show()		


# Temporal dynamics

# Demographics
#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matrices[0][0][1][:,12],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[1][0][1][:,12],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,12],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[1][0][2][:,12],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control[:,12],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.title("Number of compatible mates - demographic and genetic rescue")
plt.show()


#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matrices[0][0][1][:,2],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[1][0][1][:,2],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,2],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[1][0][2][:,2],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control[:,2],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("# vegetative individuals")
#plt.ylabel("# reproductive individuals with pollen flow")
plt.title("Number of vegetative individuals - demographic and genetic rescue")
plt.show()


#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matrices[0][0][1][:,4],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[1][0][1][:,4],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,4],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[1][0][2][:,4],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control[:,4],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("# reproductive individuals")
#plt.ylabel("# reproductive individuals with pollen flow")
plt.title("Number of reproductive individuals - demographic and genetic rescue")
plt.show()

#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matrices[0][0][1][:,6],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[1][0][1][:,6],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,6],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[1][0][2][:,6],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control[:,6],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("Age")
#plt.ylabel("# reproductive individuals with pollen flow")
plt.title("Age")
plt.show()

plt.plot(matrices[0][0][1][:,8],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[1][0][1][:,8],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,8],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[1][0][2][:,8],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control[:,8],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("Age")
#plt.ylabel("# reproductive individuals with pollen flow")
plt.title("Age")
plt.show()

#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matrices[0][0][1][:,1],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[1][0][1][:,1],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,1],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[1][0][2][:,1],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control[:,1],'r--', label = "Control")
plt.legend(loc='lower left')
plt.xlabel("Generations")
plt.ylabel("# runs persisiting to that generation")
#plt.ylabel("# reproductive individuals with pollen flow")
plt.title("# runs persisiting to that generation")
plt.show()

#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matricesfit[0][0][1][:,11],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
plt.plot(matricesfit[3][0][1][:,11],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
plt.plot(matricesfit[0][0][2][:,11],'g--', label = "pop increase = 10 - prob new S = 0.2",linewidth=2.0)
plt.plot(matricesfit[3][0][2][:,11],'k--', label = "pop increase = 40 - prob new S = 0.2", linewidth=2.0)
plt.plot(npmat_control_f[:,11],'r--', label = "Control", linewidth=2.0)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.xlabel("Generations",fontsize=16)
plt.ylabel("Average seeds per plant post dispersal",fontsize=16)
plt.title("Average seeds per plant post dispersal")
plt.show()

fig, ax = plt.subplots(1,2, figsize=(15, 6), facecolor='w', edgecolor='k')
ax[0].plot(matricesfit[0][0][0][:,11],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
ax[0].plot(matricesfit[2][0][0][:,11],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
ax[0].plot(matricesfit[0][0][1][:,11],'g--', label = "pop increase = 10 - prob new S = 0.2",linewidth=2.0)
ax[0].plot(matricesfit[2][0][1][:,11],'k--', label = "pop increase = 40 - prob new S = 0.2", linewidth=2.0)
ax[0].plot(npmat_control_f[:,11],'r--', label = "Control", linewidth=2.0)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
ax[0].legend()
ax[0].set_xlabel("Generations",fontsize=16)
ax[0].set_ylabel("Average seeds per plant post dispersal",fontsize=16)
ax[0].set_title("Average seeds per plant post dispersal")

ax[1].plot(matricesfit[0][0][0][:,11],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
ax[1].plot(matricesfit[2][0][0][:,11],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
ax[1].plot(matricesfit[0][0][5][:,11],'g--', label = "pop increase = 10 - prob new S = 0.8",linewidth=2.0)
ax[1].plot(matricesfit[2][0][5][:,11],'k--', label = "pop increase = 40 - prob new S = 0.8", linewidth=2.0)
ax[1].plot(npmat_control_f[:,11],'r--', label = "Control", linewidth=2.0)
ax[1].legend()
ax[1].set_xlabel("Generations",fontsize=16)
ax[1].set_ylabel("Average seeds per plant post dispersal",fontsize=16)
ax[1].set_title("Average seeds per plant post dispersal")
plt.show()

plt.plot(matricesfit[0][0][1][:,1],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
plt.plot(matricesfit[3][0][1][:,1],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
plt.plot(matricesfit[0][0][2][:,1],'g--', label = "pop increase = 10 - prob new S = 0.2",linewidth=2.0)
plt.plot(matricesfit[3][0][2][:,1],'k--', label = "pop increase = 40 - prob new S = 0.2", linewidth=2.0)
plt.plot(npmat_control_f[:,1],'r--', label = "Control", linewidth=2.0)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='lower left')
plt.xlabel("Generations",fontsize=16)
plt.ylabel("# pollen donors/female",fontsize=16)
plt.title("# pollen donors/female")
plt.show()

plt.plot(matricesfit[0][3][1][:,11],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesfit[3][3][1][:,11],'k', label = "pop increase = 40 - prob new S = 0")
plt.plot(matricesfit[0][3][3][:,11],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesfit[3][3][3][:,11],'k--', label = "pop increase = 40 - prob new S = 0.2")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("Average seeds per plant post dispersal")
plt.title("Average seeds per plant post dispersal - demographic and genetic rescue")
plt.show()

#plt.plot(npmat1[:,12],'r', label = "pop increase = 10")
plt.plot(matricesgen[0][0][1][:,4],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
plt.plot(matricesgen[3][0][1][:,4],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
plt.plot(matricesgen[0][0][2][:,4],'g--', label = "pop increase = 10 - prob new S = 0.2", linewidth=2.0)
plt.plot(matricesgen[3][0][2][:,4],'k--', label = "pop increase = 40 - prob new S = 0.2", linewidth=2.0)
plt.plot(npmat_control_g[:,4],'r--', label = "Control", linewidth=2.0)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='lower left')
plt.xlabel("Generations",fontsize=16)
plt.ylabel("Average # of S alleles",fontsize=16)
plt.title("Average # of S alleles ", fontsize=20)
plt.show()

fig, ax = plt.subplots(1,2, figsize=(15, 6), facecolor='w', edgecolor='k')
ax[0].plot(matricesgen[0][0][0][:,4],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
ax[0].plot(matricesgen[2][0][0][:,4],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
ax[0].plot(matricesgen[0][0][1][:,4],'g--', label = "pop increase = 10 - prob new S = 0.2",linewidth=2.0)
ax[0].plot(matricesgen[2][0][1][:,4],'k--', label = "pop increase = 40 - prob new S = 0.2", linewidth=2.0)
ax[0].plot(npmat_control_g[:,4],'r--', label = "Control", linewidth=2.0)
ax[0].tick_params(axis='y', labelsize=16)
ax[0].tick_params(axis='x', labelsize=16)
ax[0].legend(loc='lower left')
ax[0].set_xlabel("Generations",fontsize=18)
ax[0].set_ylabel("Average # of S alleles",fontsize=18)
#ax[0].set_title("Average # of S alleles")

ax[1].plot(matrices[0][0][0][:,4],'g', label = "pop increase = 10 - prob new S = 0", linewidth=2.0)
ax[1].plot(matrices[2][0][0][:,4],'k', label = "pop increase = 40 - prob new S = 0", linewidth=2.0)
ax[1].plot(matrices[0][0][1][:,4],'g--', label = "pop increase = 10 - prob new S = 0.8",linewidth=2.0)
ax[1].plot(matrices[2][0][1][:,4],'k--', label = "pop increase = 40 - prob new S = 0.8", linewidth=2.0)
ax[1].plot(npmat_control[:,4],'r--', label = "Control", linewidth=2.0)
ax[1].tick_params(axis='x', labelsize=16)
ax[1].tick_params(axis='y', labelsize=16)
ax[1].set_xlabel("Generations",fontsize=18)
ax[1].set_ylabel("# reproductive individuals",fontsize=18)
#ax[1].set_title("# reproductive individuals")
plt.show()

plt.plot(matricesgen[0][0][0][:,14],'g', label = "pop increase = 10 - prob new S = 0",linewidth=2.0)
plt.plot(matricesgen[3][0][0][:,14],'k', label = "pop increase = 40 - prob new S = 0",linewidth=2.0)
plt.plot(matricesgen[0][0][1][:,14],'g--', label = "pop increase = 10 - prob new S = 0.2",linewidth=2.0)
plt.plot(matricesgen[3][0][1][:,14],'k--', label = "pop increase = 40 - prob new S = 0.2",linewidth=2.0)
plt.plot(npmat_control_g[:,14],'r--', label = "Control")
plt.legend(loc='upper left')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel("Generations",fontsize=16)
plt.ylabel("Inbreeding coefficient (Fis)",fontsize=16)
#plt.title("Inbreeding coefficient (Fis)")
plt.show()


plt.plot(matricesgen[0][0][1][:,12],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesgen[3][0][1][:,12],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matricesgen[0][0][2][:,12],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesgen[3][0][2][:,12],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control_g[:,12],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations",fontsize=14)
plt.ylabel("He",fontsize=14)
plt.title("Average Allelic diversity ")
plt.show()


plt.plot(matricesgen[0][0][0][:,2],'g', label = "pop increase = 10 - prob new S = 0",linewidth=2.0)
plt.plot(matricesgen[3][0][0][:,2],'k', label = "pop increase = 40 - prob new S = 0",linewidth=2.0)
plt.plot(matricesgen[0][0][1][:,2],'g--', label = "pop increase = 10 - prob new S = 0.2",linewidth=2.0)
plt.plot(matricesgen[3][0][1][:,2],'k--', label = "pop increase = 40 - prob new S = 0.2",linewidth=2.0)
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,2]-matricesgen[1][0][2][:,3],matricesgen[1][0][2][:,2]+matricesgen[1][0][2][:,3])
plt.plot(npmat_control_g[:,2],'r--', label = "Control",linewidth=2.0)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc='lower left')
plt.xlabel("Generations",fontsize=16)
plt.ylabel("Average neutral alleles",fontsize=16)
#plt.title("Average neutral alleles")
#plt.title("Average neutral alleles")
plt.show()

plt.plot(matricesgen[0][0][1][:,6],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesgen[1][0][1][:,6],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matricesgen[0][0][2][:,6],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesgen[1][0][2][:,6],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control_g[:,6],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("Variance neutral alleles")
plt.title("Variance neutral alleles - demographic and genetic rescue")
plt.show()

plt.plot(matricesgen[0][0][0][:,10],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesgen[3][0][0][:,10],'k', label = "pop increase = 40 - prob new S = 0")
plt.plot(matricesgen[0][0][1][:,10],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesgen[3][0][1][:,10],'k--', label = "pop increase = 40 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
plt.plot(npmat_control_g[:,10],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
plt.legend()
plt.xlabel("Generations")
plt.ylabel("Ho")
plt.title("Ho - demographic and genetic rescue")
plt.show()

## Spatial rescue
plt.plot(matrices[0][0][1][:,10],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matrices[0][2][1][:,10],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matrices[0][0][2][:,10],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matrices[0][2][2][:,10],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.plot(npmat_control_g[:,10],'r--', label = "Control")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("Ho")
plt.title("Ho - demographic and genetic rescue")
plt.show()

plt.plot(matricesgen[0][0][1][:,4],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesgen[0][2][1][:,4],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matricesgen[0][0][2][:,4],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesgen[0][2][2][:,4],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("# Pollen donors / female")
plt.title("Number of compatible mates - demographic and genetic rescue")
plt.show()

plt.plot(matricesfit[0][0][1][:,7],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesfit[0][2][1][:,7],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matricesfit[0][0][2][:,7],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesfit[0][2][2][:,7],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("# Pollen donors / female")
plt.title("Number of compatible mates - demographic and genetic rescue")
plt.show()

plt.plot(matricesfit[0][0][1][:,13],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesfit[0][2][1][:,13],'k', label = "pop increase = 20 - prob new S = 0")
plt.plot(matricesfit[0][0][2][:,13],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesfit[0][2][2][:,13],'k--', label = "pop increase = 20 - prob new S = 0.2")
plt.legend()
plt.xlabel("Generations")
plt.ylabel("# Pollen donors / female")
plt.title("Number of compatible mates - demographic and genetic rescue")
plt.show()



-----------------------------------------------------------------------

popsize = ['10','30']
safesites = ['0.0','0.35']
probnews = ['0.0','0.2']


# demographic files 

persistn = np.zeros((len(popsize),len(safesites),len(probnews)))
counti = 0
matrices_neutral = []        
for i in popsize:
	matj = []
	countj = 0
	for j in safesites:
		mati = []
		countk = 0
		for k in probnews:
			mat1 = []
			file1 = open('demo_neutral_' + i + j + k + '.dat', 'rb')
			reader1 = csv.reader(file1, delimiter='\t')
			for row in reader1:
				mat1.append(row)
			persistence = mat1[0:11]
			if persistence:
				persistn[counti,countj,countk] = float(persistence[7][0].strip("Average Persistence: "))	
			del mat1[0:11]
			mat11 =[]
			for l in mat1:
				if len(l) > 1:
					mat11.append(l)        
			npmat1 = np.zeros((len(mat11),14))
			for ii in range(len(mat11)):
				for jj in range(14):
					if mat11[ii][jj] == '.':
						mat11[ii][jj] = 0
					npmat1[ii,jj] = float(mat11[ii][jj])
			file1.close() 
			mati.append(npmat1)
			countk += 1
		countj += 1
		matj.append(mati)
	matrices_neutral.append(matj)
	counti += 1

# fitness files

counti = 0
matricesfitn = []        
for i in popsize:
	matj = []
	countj = 0
	for j in safesites:
		mati = []
		countk = 0
		for k in probnews:
			mat1 = []
			file1 = open('fit_neutral_' + i + j + k + '.dat', 'rb')
			reader1 = csv.reader(file1, delimiter='\t')
			for row in reader1:
				mat1.append(row)
			del mat1[0:11]
			mat11 =[]
			for l in mat1:
				if len(l) > 1:
					mat11.append(l)        
			npmat1 = np.zeros((len(mat11),17))
			for ii in range(len(mat11)):
				for jj in range(17):
					if mat11[ii][jj] == '.':
						mat11[ii][jj] = 0
					npmat1[ii,jj] = float(mat11[ii][jj])
			file1.close() 
			mati.append(npmat1)
			countk += 1
		countj += 1
		matj.append(mati)
	matricesfitn.append(matj)
	counti += 1


# genetics files

counti = 0
matricesgenn = []        
for i in popsize:
	matj = []
	countj = 0
	for j in safesites:
		mati = []
		countk = 0
		for k in probnews:
			mat1 = []
			file1 = open('gen_neutral_' + i + j + k + '.dat', 'rb')
			reader1 = csv.reader(file1, delimiter='\t')
			for row in reader1:
				mat1.append(row)
			del mat1[0:8]
			mat11 =[]
			for l in mat1:
				if len(l) > 1:
					mat11.append(l)        
			npmat1 = np.zeros((len(mat11),16))
			for ii in range(len(mat11)):
				for jj in range(16):
					if mat11[ii][jj] == '.':
						mat11[ii][jj] = 0
					npmat1[ii,jj] = float(mat11[ii][jj])
			file1.close() 
			mati.append(npmat1)
			countk += 1
		countj += 1
		matj.append(mati)
	matricesgenn.append(matj)
	counti += 1


## Neutral alleles prob increase simulations (without spatial rescue ->> site increase = 0.0) 
plt.plot(matricesgenn[0][0][0][:,12],'g', label = "pop increase = 10 - prob new S = 0",linewidth=2.0)
plt.plot(matricesgenn[1][0][0][:,12],'k', label = "pop increase = 30 - prob new S = 0")
plt.plot(matricesgenn[0][0][1][:,12],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesgenn[1][0][1][:,12],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
plt.plot(npmat_control_g[:,12],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
plt.legend()
plt.xlabel("Generations")
plt.ylabel("He")
plt.title("He - demographic and genetic rescue - Site inc = 0.0")
plt.show()

plt.plot(matricesgenn[0][1][0][:,12],'g', label = "pop increase = 10 - prob new S = 0")
plt.plot(matricesgenn[1][1][0][:,12],'k', label = "pop increase = 30 - prob new S = 0")
plt.plot(matricesgenn[0][1][1][:,12],'g--', label = "pop increase = 10 - prob new S = 0.2")
plt.plot(matricesgenn[1][1][1][:,12],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
plt.plot(npmat_control_g[:,12],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
plt.legend()
plt.xlabel("Generations")
plt.ylabel("He")
plt.title("He - demographic and genetic rescue - Site inc = 0.35")
plt.show()

## Neutral alleles prob increase simulations (without spatial rescue ->> site increase = 0.0)
fig, axs = plt.subplots(1,2, figsize=(15, 6), facecolor='w', edgecolor='k') 
axs[0].plot(matricesgenn[0][0][0][:,14],'g', label = "pop increase = 10 - prob new S = 0",linewidth=2.0)
axs[0].plot(matricesgenn[1][0][0][:,14],'k', label = "pop increase = 30 - prob new S = 0")
axs[0].plot(matricesgenn[0][0][1][:,14],'g--', label = "pop increase = 10 - prob new S = 0.2")
axs[0].plot(matricesgenn[1][0][1][:,14],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
axs[0].plot(npmat_control_g[:,14],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
axs[0].legend()
axs[0].set_xlabel("Generations")
axs[0].set_ylabel("Fis")
axs[0].set_title("Fis - demographic and genetic rescue - Site inc = 0.0")

axs[1].plot(matricesgenn[0][1][0][:,14],'g', label = "pop increase = 10 - prob new S = 0")
axs[1].plot(matricesgenn[1][1][0][:,14],'k', label = "pop increase = 30 - prob new S = 0")
axs[1].plot(matricesgenn[0][1][1][:,14],'g--', label = "pop increase = 10 - prob new S = 0.2")
axs[1].plot(matricesgenn[1][1][1][:,14],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
axs[1].plot(npmat_control_g[:,14],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
axs[1].legend()
axs[1].set_xlabel("Generations")
axs[1].set_ylabel("Fis")
axs[1].set_title("Fis - demographic and genetic rescue - Site inc = 0.35")
plt.show()


## Neutral alleles prob increase simulations (without spatial rescue ->> site increase = 0.0)
fig, axs = plt.subplots(1,2, figsize=(15, 6), facecolor='w', edgecolor='k') 
axs[0].plot(matricesgenn[0][0][0][:,10],'g', label = "pop increase = 10 - prob new S = 0",linewidth=2.0)
axs[0].plot(matricesgenn[1][0][0][:,10],'k', label = "pop increase = 30 - prob new S = 0")
axs[0].plot(matricesgenn[0][0][1][:,10],'g--', label = "pop increase = 10 - prob new S = 0.2")
axs[0].plot(matricesgenn[1][0][1][:,10],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
axs[0].plot(npmat_control_g[:,10],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
axs[0].legend()
axs[0].set_xlabel("Generations")
axs[0].set_ylabel("Ho")
axs[0].set_title("Ho - demographic and genetic rescue - Site inc = 0.0")

axs[1].plot(matricesgenn[0][1][0][:,10],'g', label = "pop increase = 10 - prob new S = 0")
axs[1].plot(matricesgenn[1][1][0][:,10],'k', label = "pop increase = 30 - prob new S = 0")
axs[1].plot(matricesgenn[0][1][1][:,10],'g--', label = "pop increase = 10 - prob new S = 0.2")
axs[1].plot(matricesgenn[1][1][1][:,10],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
axs[1].plot(npmat_control_g[:,10],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
axs[1].legend()
axs[1].set_xlabel("Generations")
axs[1].set_ylabel("Ho")
axs[1].set_title("Ho - demographic and genetic rescue - Site inc = 0.35")
plt.show()

## Neutral alleles prob increase simulations (without spatial rescue ->> site increase = 0.0)
fig, axs = plt.subplots(1,2, figsize=(15, 6), facecolor='w', edgecolor='k') 
axs[0].plot(matricesgenn[0][0][0][:,14],'g', label = "pop increase = 10 - prob new S = 0",linewidth=2.0)
axs[0].plot(matricesgenn[1][0][0][:,14],'k', label = "pop increase = 30 - prob new S = 0")
axs[0].plot(matricesgenn[0][0][1][:,14],'g--', label = "pop increase = 10 - prob new S = 0.2")
axs[0].plot(matricesgenn[1][0][1][:,14],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
axs[0].plot(npmat_control_g[:,14],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
axs[0].legend()
axs[0].set_xlabel("Generations")
axs[0].set_ylabel("Fis")
axs[0].set_title("Fis - demographic and genetic rescue - Site inc = 0.0")

axs[1].plot(matricesgenn[0][0][0][:,10],'g', label = "pop increase = 10 - prob new S = 0")
axs[1].plot(matricesgenn[1][0][0][:,10],'k', label = "pop increase = 30 - prob new S = 0")
axs[1].plot(matricesgenn[0][0][1][:,10],'g--', label = "pop increase = 10 - prob new S = 0.2")
axs[1].plot(matricesgenn[1][0][1][:,10],'k--', label = "pop increase = 30 - prob new S = 0.2")
#plt.fill_between(range(0,500),matricesgen[1][0][2][:,10]-matricesgen[1][0][2][:,11],matricesgen[1][0][2][:,10]+matricesgen[1][0][2][:,11])
axs[1].plot(npmat_control_g[:,10],'r--', label = "Control")
#plt.fill_between(range(0,500), npmat_control_g[:,10]-npmat_control_g[:,11], npmat_control_g[:,10]+npmat_control_g[:,11])
axs[1].legend()
axs[1].set_xlabel("Generations")
axs[1].set_ylabel("Fis")
axs[1].set_title("Fis - demographic and genetic rescue - Site inc = 0.35")
plt.show()
-----------------------------------------------------------------------

import numpy as np
import csv
import matplotlib.pyplot as plt

fileg1 = open("demo_rescue_10_gen.dat",'rb')
readerg1 = csv.reader(fileg1, delimiter='\t')

fileg2 = open("spatial_rescue_02_gen.dat",'rb')
readerg2 = csv.reader(fileg2, delimiter='\t')

fileg3 = open("spatial_rescue_03_gen.dat",'rb')
readerg3 = csv.reader(fileg3, delimiter='\t')

matgen1 = []
for row in readerg1:
    matgen1.append(row)

matgen2 = []
for row in readerg2:
    matgen2.append(row)

matgen3 = []
for row in readerg3:
    matgen3.append(row)

fileg1.close()
fileg2.close()
fileg3.close()

del matgen1[0:7]    
matgen11 = []

for i in matgen1:
    if len(i) > 1:
        matgen11.append(i)

del matgen2[0:7]    
matgen22 = []

for i in matgen2:
    if len(i) > 1:
        matgen22.append(i)

del matgen3[0:7]    
matgen33 = []

for i in matgen3:
    if len(i) > 1:
        matgen33.append(i)

npmatg1 = np.zeros((len(matgen11),len(matgen11[0])))
npmatg2 = np.zeros((len(matgen22),len(matgen11[0])))
npmatg3 = np.zeros((len(matgen33),len(matgen11[0])))
        
for i in range(len(mat11)):
   for j in range(len(matgen11[0])):
       npmatg1[i,j] = float(matgen11[i][j])
       npmatg2[i,j] = float(matgen22[i][j])
       npmatg3[i,j] = float(matgen33[i][j])


#plt.plot(npmatg1[:,10],'r', label = "pop increase = 10")
plt.plot(npmatg2[:,10],'g', label = "spatial increase = 0.2")
plt.plot(npmatg3[:,10],'k', label = "spatial increase = 0.3")
plt.legend(loc = "lower left")
plt.xlabel("Generations")
plt.ylim(ymax = 0.7, ymin = 0.5)
plt.ylabel("Observed heterozygosity (Ho) - spatial rescue")
plt.show()


#plt.plot(npmatg1[:,14],'r', label = "pop increase = 10")
plt.plot(npmatg2[:,14],'g', label = "spatial increase = 0.2")
plt.plot(npmatg3[:,14],'k', label = "spatial increase = 0.3")
plt.legend(loc = "lower left")
plt.xlabel("Generations")
#plt.ylim(ymax = 0.7, ymin = 0.5)
plt.ylabel("Inbreeding (Fis) - spatial rescue")
plt.show()


-------------------------------------------------------------



import numpy as np
import csv
import matplotlib.pyplot as plt

filef1 = open("demo_rescue_10_fit.dat",'rb')
readerf1 = csv.reader(filef1, delimiter='\t')

filef2 = open("spatial_rescue_02_fit.dat",'rb')
readerf2 = csv.reader(filef2, delimiter='\t')

filef3 = open("spatial_rescue_03_fit.dat",'rb')
readerf3 = csv.reader(filef3, delimiter='\t')

matfit1 = []
for row in readerf1:
    matfit1.append(row)

matfit2 = []
for row in readerf2:
    matfit2.append(row)

matfit3 = []
for row in readerf3:
    matfit3.append(row)

filef1.close()
filef2.close()
filef3.close()

del matfit1[0:8]    
matfit11 = []

for i in matfit1:
    if len(i) > 1:
        matfit11.append(i)

del matfit2[0:8]    
matfit22 = []

for i in matfit2:
    if len(i) > 1:
        matfit22.append(i)

del matfit3[0:8]    
matfit33 = []

for i in matfit3:
    if len(i) > 1:
        matfit33.append(i)

npmatf1 = np.zeros((len(matfit11),len(matfit11[0])))
npmatf2 = np.zeros((len(matfit22),len(matfit11[0])))
npmatf3 = np.zeros((len(matfit33),len(matfit11[0])))
        
for i in range(len(matfit11)):
   for j in range(len(matfit11[0])):
       npmatf1[i,j] = float(matfit11[i][j])
       npmatf2[i,j] = float(matfit22[i][j])
       npmatf3[i,j] = float(matfit33[i][j])


#plt.plot(npmatf1[:,1],'r', label = "pop increase = 10")
plt.plot(npmatf2[:,1],'g', label = "spatial increase = 0.2")
plt.plot(npmatf3[:,1],'k', label = "spatial increase = 0.3")
plt.legend(loc = "lower left")
plt.xlabel("Generations")
#plt.ylim(ymax = 0.7, ymin = 0.5)
plt.ylabel("Number of pollen donors - spatial rescue")
plt.show()


#plt.plot(npmatf1[:,5],'r', label = "pop increase = 10")
plt.plot(npmatf2[:,5],'g', label = "spatial increase = 0.2")
plt.plot(npmatf3[:,5],'k', label = "spatial increase = 0.3")
plt.legend(loc = "lower left")
plt.xlabel("Generations")
#plt.ylim(ymax = 0.7, ymin = 0.5)
plt.ylabel("Average seeds before dispersal - spatial rescue")
plt.show()


#plt.plot(npmatf1[:,11],'r', label = "pop increase = 10")
plt.plot(npmatf2[:,11],'g', label = "spatial increase = 0.2")
plt.plot(npmatf3[:,11],'k', label = "spatial increase = 0.3")
plt.legend()
plt.xlabel("Generations")
#plt.ylim(ymax = 0.7, ymin = 0.5)
plt.ylabel("Average seeds after dispersal - spatial rescue")
plt.show()
