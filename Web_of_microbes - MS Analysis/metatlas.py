from matplotlib import pyplot as plt
import requests, json
import numpy as np
import os
import metatlas
from scipy.optimize import leastsq
from math import exp
import platform
import csv
import time
import re
import sys
import copy
from bokeh.plotting import *

def export_peakData_to_spreadsheet(filename,export_fileIds,fileInfo,data,dictData):
	import csv
	export_filenames = []
	for i,myFile in enumerate(export_fileIds):
	    for j,fid in enumerate(fileInfo['fid']):
	        if fid == myFile:
	            export_filenames.append(fileInfo['name'][j])
	fid = open(filename,'wb')
	fid.write('%s\t' % 'compound')
	for filename in export_filenames:
	    fid.write('%s\t' % filename)
	fid.write('\n')
	for i,datum in enumerate(data):
	    fid.write('%s\t' % dictData[u'compounds'][i]['name'])
	    mz = float(dictData[u'compounds'][i][u'mz'])
	    mzTol = float(dictData[u'compounds'][i][u'mz_threshold'])
	    mzMin = mz - mz*mzTol/1.0e6
	    mzMax = mz + mz*mzTol/1.0e6
	    rtMin = float(dictData[u'compounds'][i][u'rt_min'])
	    rtMax = float(dictData[u'compounds'][i][u'rt_max'])
	    for j,myFile in enumerate(export_fileIds):
	        idx = np.logical_and( datum[:,2]==myFile, datum[:,0]>=rtMin, datum[:,0]<=rtMax )
	        if np.sum(idx)>0:
	            x1 = datum[:,0][idx]
	            y1 = datum[:,1][idx]
	            # y1 = y1 - np.min(y1)
	            myname = dictData[u'compounds'][i]['name']
	            if myname.startswith('ist'):
	                y1 = y1[:]
	            else:    
	                y1 = y1[:] / fileInfo['normalization_factor'][j]
	            fid.write('%5.2f\t' % np.sum(y1))
	        else:
	            fid.write('%5.2f\t' % 0)
	    fid.write('\n')
	fid.close()

def shareExperiments(allUsers,allPerms,client, myExperimentID):
    for aUser in allUsers:
        payload = {"user":aUser,"perms":allPerms}
        sendData=json.dumps(payload)
        # print sendData
        url = 'https://metatlas.nersc.gov/api/experiment/%s/share/' % myExperimentID
        # print url
        r = client.post(url, data=sendData)
        return r
        # print r.content
		
def listMyExperiments(client):
	url = 'https://metatlas.nersc.gov/api/experiment'
	r = client.get(url)
	experiments = json.loads(r.content)
	for experiment in experiments:
		print 'NAME = %s;ID = %s' % (experiment[u'name'], experiment[u'_id'])
	return experiments

def exportfilelist(myExperimentID, client, wd, filename):
	url = 'https://metatlas.nersc.gov/api/experiment/%s' % myExperimentID
	r = client.get(url)
	files = json.loads(r.content)
	fileInfo = {'fid':[],'name':[],'status':[]};
	print "There are %s files in this experiment. Open the file: \n%s\\%s \nAssign groups, polarity, plot order and \
normalization factors; save as a txt file." % (len(files['runs']), wd[0],filename)
	fid = open(filename,'wb')
	fid.write('index\tstatus\tname\tfid\tpolarity\tgroup\tinclusion_order\tnormalization_factor\n')
	for i,myRun in enumerate(files[u'runs']):
		splitPathToFile = os.path.split(myRun[u'in_file'])
		fid.write('%d\t%d\t%s\t%d\tpos\tgroup1\n' % (i,myRun[u'pending'],splitPathToFile[1],myRun[u'_id'][u'file_id']))
		if myRun[u'pending'] == 0:
			fileInfo['fid'].append(myRun[u'_id'][u'file_id'])
			fileInfo['name'].append(splitPathToFile[1])
			fileInfo['status'].append(myRun[u'pending']) #only keep if status is 0
	pathYouWant = splitPathToFile[0]
	return files[u'runs'][0][u'_id'][u'array_name']
	fid.close()
	
def getExpName(myExperimentID, experiments):
	for experiment in experiments:
		if myExperimentID == experiment[u'_id']:
			myExpName= experiment[u'name']
	print "Experiment loaded: %s" % (myExpName)
	return myExpName
	
def uploadgroupings(filename):
	import csv
	with open(filename,'rU') as file_object:
		newfileInfo = list(csv.DictReader(file_object, dialect='excel-tab'))
	keys = newfileInfo[0].iterkeys()
	fileInfo = {key: [d[key] for d in newfileInfo] for key in keys}
	fileInfo['fid'] = map(int, fileInfo['fid'])
	fileInfo['index'] = map(int, fileInfo['index'])
	fileInfo['inclusion_order'] = map(int, fileInfo['inclusion_order'])
	fileInfo['status'] = map(int, fileInfo['status'])
	fileInfo['normalization_factor'] = map(float, fileInfo['normalization_factor'])
	print "%s - first file" % fileInfo['name'][:1]
	return fileInfo
	
def getfileIDs(fileInfo,polarity):
	idx = np.argsort(fileInfo['inclusion_order'])
	#print idx
	export_fileIds = np.asarray(fileInfo['fid'])[idx]
	# print export_fileIds
	# print fileInfo['inclusion_order']
	# print idx
	print "%s total files IDs" % len(export_fileIds)
	if polarity == 1:
		export_fileIds_pos = []
		idx_pos=[]
		for j,each in enumerate(fileInfo['polarity']):
			if each == "pos":
				export_fileIds_pos.append(fileInfo['fid'][j])
				idx_pos.append(fileInfo['inclusion_order'][j])
		idx_pos_sort=np.argsort(idx_pos)
		export_fileIds_pos=np.asarray(export_fileIds_pos)[idx_pos_sort]
		print "%s POS files IDs" % len(export_fileIds_pos)
		return export_fileIds_pos
	elif polarity == 0:
		export_fileIds_neg = []
		idx_neg=[]
		for j,each in enumerate(fileInfo['polarity']):
			if each == "Neg":
				export_fileIds_neg.append(fileInfo['fid'][j])
				idx_neg.append(fileInfo['inclusion_order'][j])
		idx_neg_sort=np.argsort(idx_neg)
		export_fileIds_neg=np.asarray(export_fileIds_neg)[idx_neg_sort]
		print "%s NEG files IDs" % len(export_fileIds_neg)
		return export_fileIds_neg
	else:
		print "error"
		
def polaritycheck(selectedpolarity):
	polarity=selectedpolarity.lower()
	if polarity == "negative":
		polarity=0
	if polarity == "positive":
		polarity=1
	print "polarity set"
	return polarity

def emptyatlas(filename):
	ColHeadings = ['name','pubchem_id','formula','neutral_mass','mz','mz_threshold','adducts','rt_max','rt_min','rt_peak']
	fid = open(filename,'wb')
	for each in ColHeadings:
		fid.write('%s\t' % each)
	fid.write('\n')
	fid.close()
	print "Empty atlas template with column labels (%s) saved to working directory" % (filename)
		
def getAtlasList(client):
	url = 'https://metatlas.nersc.gov/api/dict/'
	r = client.get(url)
	allAtlases = json.loads(r.text)
	for atlas in allAtlases:
		atlas_str = '%s has an atlas named %s has the ID: %s' % (atlas[u'creator'],atlas['name'], atlas['_id'])
		print atlas_str
	return allAtlases
	
def getAtlasEntries(atlasID, client, atlasName):
	url = 'https://metatlas.nersc.gov/api/dict/%s/' % atlasID
	r = client.get(url)
	dictData = json.loads(r.text)
	print "%s entries in atlas: %s" % (len(dictData[u'compounds']), atlasName)
	return dictData

def getAtlasName(dictID, allAtlases):

	for atlas in allAtlases:
		if dictID == atlas['_id']:
			atlasName=atlas['name']
	print "Selected atlas: %s" % (atlasName)
	return atlasName
	
def addCompAtlas(filename, client, dictId, atlasName):
	with open(filename,'rU') as file_object:
		payload = list(csv.DictReader(file_object, dialect='excel-tab'))
	url = 'https://metatlas.nersc.gov/api/dict/%s/' % dictId
	r = client.post(url, data=json.dumps(payload))
	print "%s compounds added to atlas: %s" %(len(payload), atlasName)
	
def exportAtlas(filename, atlasName, dictData):
	myList = ['name','pubchem_id','formula','neutral_mass','mz','mz_threshold','adducts','rt_max','rt_min','rt_peak']
	fid = open(filename,'wb')
	for listItem in myList:
		fid.write('%s\t' % listItem)
	fid.write('\n')
	for i,compound in enumerate(dictData[u'compounds']):
		for listItem in myList:
			fid.write('%s\t' % compound[listItem])
		fid.write('\n')
	fid.close()
	print "%s compounds exported from atlas %s (%s)" % (len(dictData[u'compounds']),atlasName,filename)
	
def getEICdata(polarity, expandRTrangeby, comp_start, comp_end, dictData, export_fileIds, exportedfiles, client):
	if comp_start == None:
		comp_start = 0
	data = []
	if comp_start==None:
		i=0
	else:
		i = comp_start
	for x in range(0,comp_start):
		data.append([])
	for compound in dictData[u'compounds'][comp_start:comp_end]:
		try: 
			extraTime=0.0
			RTrange=float(compound[u'rt_max'])-float(compound[u'rt_min'])
			if RTrange < 1.0:
				expandto1min =(1.0-RTrange)/2.0
			else:
				expandto1min = 0.0
			extraTime=expandRTrangeby+expandto1min
			newdata=metatlas.getEICForCompounds(compound,exportedfiles,export_fileIds,extraTime,client,polarity)
			data.append(newdata)
			if len(newdata)>0:
				print i, " EIC obtained for ", compound[u'name'], "RTrange: %s, extraTime: %s, Expandedrange: %s" % (RTrange, extraTime, RTrange+extraTime*2)
			else:
				print i, " no data for EIC for ", compound[u'name'], "RTrange: %s, extraTime: %s, Expandedrange: %s" % (RTrange, extraTime, RTrange+extraTime*2)
		except:
			data.append([])
			print i, " Error obtain EIC for ", compound[u'name']
		i = i+1
		time.sleep(0.5)
	if len(dictData[u'compounds'])+1> comp_end:
		for x in range(comp_end-1,len(dictData[u'compounds'])+1):
			data.append([])
	return data
	
def authenticateUser(userFile):
    authURL = 'https://metatlas.nersc.gov/client/login/'
    file = open(userFile, 'r')
    userID = file.readline()[:-1]
    if platform.system()=="Windows":
        userPassword = file.readline()[:]
    else:
        userPassword = file.readline()[:-1]
    file.close()

    client = requests.Session()
    # Retrieve the CSRF token first
    client.get(authURL)  # sets cookie
    csrftoken = client.cookies['csrftoken']
    login_data = dict(username=userID, password=userPassword, csrfmiddlewaretoken=csrftoken, next='/')
    r = client.post(authURL, data=login_data, headers=dict(Referer=authURL))
    return client

def getEICForCompounds_oneLighter(compound,myArray,files_I_want,rtTol,client,polarity):
	if isinstance(files_I_want,int):
		myList = str(files_I_want)
	else:
		myList = ','.join(map(str, files_I_want))
	mz = float(compound[u'mz'])
	mzTol = float(compound[u'mz_threshold'])
	mzMin = mz - mz*mzTol/1.0e6 - 1.003355
	mzMax = mz + mz*mzTol/1.0e6 - 1.003355
	rtMin = float(compound[u'rt_min'])-rtTol
	rtMax = float(compound[u'rt_max'])+rtTol
	rtPeak = float(compound[u'rt_peak'])

	payload = {'L':1,'P':polarity,'arrayname':myArray,'fileidlist':myList,
	          'max_mz':mzMax,'min_mz':mzMin,
	          'min_rt':rtMin,'max_rt':rtMax,
	          'nsteps':20000,'queryType':'XICofFile_mf'}
	url = 'https://metatlas.nersc.gov/api/run'
	r = client.get(url,params=payload)
	data = np.asarray(json.loads(r.content))
	return data

def getEICForCompounds(compound,myArray,files_I_want,rtTol,client,polarity):
	if isinstance(files_I_want,int):
		myList = str(files_I_want)
	else:
		myList = ','.join(map(str, files_I_want))
	mz = float(compound[u'mz'])
	mzTol = float(compound[u'mz_threshold'])
	mzMin = mz - mz*mzTol/1.0e6
	mzMax = mz + mz*mzTol/1.0e6
	rtMin = float(compound[u'rt_min'])-rtTol
	rtMax = float(compound[u'rt_max'])+rtTol
	rtPeak = float(compound[u'rt_peak'])

	payload = {'L':1,'P':polarity,'arrayname':myArray,'fileidlist':myList,
	          'max_mz':mzMax,'min_mz':mzMin,
	          'min_rt':rtMin,'max_rt':rtMax,
	          'nsteps':20000,'queryType':'XICofFile_mf'}
	url = 'https://metatlas.nersc.gov/api/run'
	r = client.get(url,params=payload)
	data = np.asarray(json.loads(r.content))
	return data
	# for myFile in files_I_want:
	#     x1 = data[:,0][(data[:,2]==myFile)]
	#     y1 = data[:,1][(data[:,2]==myFile)]
	#     idx = np.argsort(x1)
	#     plt.plot(x1[idx],y1[idx])
	# plt.xlabel('Time (min)')
	# plt.ylabel('TIC Intensity (au)')

def getEICForCompound(compound,myArray,runId,rtTol,client,polarity):
    #polarity is 1 for pos and 0 for neg
	mz = float(compound[u'mz'])
	mzTol = float(compound[u'mz_threshold'])
	mzMin = mz - mz*mzTol/1.0e6
	mzMax = mz + mz*mzTol/1.0e6
	rtMin = float(compound[u'rt_min'])
	rtMax = float(compound[u'rt_max'])
	rtPeak = float(compound[u'rt_peak'])
	payload = {'L':1,'P':polarity,'arrayname':myArray,'fileid':runId,
	'max_mz':mzMax,'min_mz':mzMin,
	'nsteps':10000,'queryType':'XICofFile'}
	url = 'https://metatlas.nersc.gov/api/run'
	r = client.get(url,params=payload)
	data = np.asarray(json.loads(r.content))
	xdata    = data[abs(data[:,0]-rtPeak)<rtTol,0]
	ydata    = data[abs(data[:,0]-rtPeak)<rtTol,1]
	peakArea = data[(data[:,0]>rtMin) & (data[:,0]<rtMax),1]
	if len(peakArea)>0:
		peakArea = sum(peakArea)
	else:
		peakArea = 0
	if len(xdata)>0:
		iMin = min(ydata)
		ydata = ydata - iMin
		iMax = max(ydata) + iMin
		ydata = ydata / iMax
	else:
		iMax = 0
	return {'eic':data,'xdata':xdata,'ydata':ydata,'name':compound[u'name'],'iMax':iMax,'peakArea':peakArea}

def createChromatogramPlots(data,compound,fitResult,ax):
	ax.plot(data['xdata'],data['ydata']*data['iMax'],'k-',data['xdata'], fitfunc(fitResult, data['xdata'])*data['iMax'],'r-',linewidth=2.0)
	ax.axvline(float(compound[u'rt_min']),linewidth=2, color='k') #original rtMin
	ax.axvline(float(compound[u'rt_max']),linewidth=2, color='k') #original rtMax
	ax.axvline(float(compound[u'rt_peak']),linewidth=2, color='g',alpha=0.5) #original rtPeak
	#     ax.axvline(x=compound[u'rt_peak'],linewidth=2, color='b') #original rtPeak
	ax.axvline(x=fitResult[1],linewidth=2, color='r') #new rtPeak
	ax.axvspan(fitResult[1]-fitResult[3]*2, fitResult[1]+fitResult[2]*2, facecolor='c', alpha=0.5) #new rtBounds
	ax.set_xlabel('Time (min)')
	ax.set_ylabel('Intensity (au)')
	ax.set_title(compound[u'name'])

def createChromatogramPlots_dataOnly(data,compound,ax):
	ax.plot(data['xdata'],data['ydata']*data['iMax'],'k-',linewidth=2.0)
	ax.axvline(float(compound[u'rt_min']),linewidth=2, color='r',alpha=0.5) #original rtMin
	ax.axvline(float(compound[u'rt_max']),linewidth=2, color='r',alpha=0.5) #original rtMax
	ax.axvline(float(compound[u'rt_peak']),linewidth=2, color='g',alpha=0.5) #original rtPeak
	#     ax.axvline(x=compound[u'rt_peak'],linewidth=2, color='b') #original rtPeak
	# ax.axvline(x=fitResult[1],linewidth=2, color='r') #new rtPeak
	# ax.axvspan(fitResult[1]-fitResult[3]*2, fitResult[1]+fitResult[2]*2, facecolor='c', alpha=0.5) #new rtBounds
	ax.set_xlabel('Time (min)')
	ax.set_ylabel('Intensity (au)')
	ax.set_title(compound[u'name'])

def fitACompound(compound,data):
	rtPeak = float(compound[u'rt_peak'])
	rtMin = float(compound[u'rt_min'])
	rtMax = float(compound[u'rt_max'])
	init  = [1.0, rtPeak, 0.1,0.1]
	out   = leastsq( errfunc, init, args=(data['xdata'], data['ydata'], rtPeak, rtMin, rtMax))
	fitResult = out[0]
	fitResult[2] = abs(fitResult[2])
	fitResult[3] = abs(fitResult[3])
	return fitResult

def fitfunc(p,x):
	a = np.zeros(x.shape)
	idx1 = x>=p[1]
	if len(idx1)>0:
		a[idx1]=p[0]*np.exp(-0.5*((x[idx1]-p[1])/p[2])**2)
	idx2 = x<p[1]
	if len(idx2)>0:
		a[idx2]=p[0]*np.exp(-0.5*((x[idx2]-p[1])/p[3])**2)
	return a

def errfunc(p,x,y,rtPeak, rtMin, rtMax):
	if (abs(p[2]) > 1.5) or (abs(p[3]) > 1.5) or (abs(p[2]) < 0.001) or (abs(p[3]) < 0.001) or (p[1] > rtMax) or (p[1] < rtMin):
		return 1e100
	else:
		# return (y-fitfunc(p,x))**2
		# idx = x > rtMin and x < rtMax
		return np.multiply((y-fitfunc(p,x))**2,np.exp(-0.5*((x-rtPeak)/0.2)**2))

def exportfilenames(dictData, filenameprefix, comp_start, comp_end):
	export_filenames = []
	for i,compound in enumerate(dictData[u'compounds'][comp_start:comp_end]):
		export_filenames.append('%s%s%s' % (filenameprefix,re.sub('[^A-Za-z0-9]+', '', str(compound[u'name'])),'.png'))
	return export_filenames
	
def subplottitles(export_fileIds,fileInfo):
	subplot_titles = []
	for i,myFile in enumerate(export_fileIds):
		for j,fid in enumerate(fileInfo['fid']):
			if fid == myFile:
				subplot_titles.append(fileInfo['name'][j].replace('.mzML',''))
	return subplot_titles
		
def makingimagefilesofplots(export_fileIds,dictData,filenameprefix, data, fileInfo,rules,polarity, comp_start, comp_end):
	export_filenames =metatlas.exportfilenames(dictData, filenameprefix,comp_start, comp_end)
	#print export_filenames
	subplot_titles = metatlas.subplottitles(export_fileIds,fileInfo)
	numCols = 12.
	nRows = int(np.ceil(len(export_fileIds)/numCols))
	if comp_start == None:
		comp_start = 0
	startedit=0 #t=len(dictData[u'compounds'][comp_start:comp_end])-len(data)
	for i,compound in enumerate(dictData[u'compounds'][comp_start:comp_end]):
		if len(data[i])>0:
			fig, ax = plt.subplots(nRows, int(numCols),figsize=(5*numCols,nRows * 5))
		#     fig.subplots_adjust(bottom=0.0, left=0.0)
			min_x_val = 1000000
			max_x_val = 0
			for j,a in enumerate(ax.flat):
				if j<len(export_fileIds):
					a.set_title(subplot_titles[j])
					x1 = data[i][:,0][(data[i][:,2]==export_fileIds[j])]
					y1 = data[i][:,1][(data[i][:,2]==export_fileIds[j])]
					y1 = y1[:] / fileInfo['normalization_factor'][j]
					idx = np.argsort(x1)
					x1 = x1[idx]
					y1 = y1[idx]
					if len(x1)>10:
						m = np.max(y1)
						y1 = y1 / m
						tempData = {'xdata':x1,'ydata':y1,'name':subplot_titles,'iMax':m}
						tempComp = copy.deepcopy(compound)
						lowercomp="".join([x.lower() for x in tempComp[u'name']])
						Timea=timeadjustments(tempComp,j,lowercomp,rules,polarity)
						tempComp[u'rt_peak']=float(tempComp[u'rt_peak'])+float(Timea)
						tempComp[u'rt_min']=float(tempComp[u'rt_min'])+float(Timea)
						tempComp[u'rt_max']=float(tempComp[u'rt_max'])+float(Timea)
						fitResult = metatlas.fitACompound(tempComp,tempData)
						metatlas.createChromatogramPlots(tempData,tempComp,fitResult,a)
						a.set_title(subplot_titles[j])
						a.set_ylim([0,np.max(data[i][:,1])])
						if np.min(data[i][:,0])<min_x_val:
							min_x_val = np.min(data[i][:,0])
						if np.max(data[i][:,0])>max_x_val:
							max_x_val = np.max(data[i][:,0])
			for j,a in enumerate(ax.flat):
				a.set_xlim([min_x_val,max_x_val])
			fig.tight_layout()        
			fig.savefig(export_filenames[i])
			fig.clear()
			print int(i)+int(comp_start), " png created for ", compound[u'name']
		else:
			print int(i)+int(comp_start), " no data for png for ",compound[u'name']

def timeadjustments(tempComp,j,lowercomp,rules,polarity):
	exec rules
	return Timeadjust
		
def outputofEICdatatotxt(export_fileIds, fileInfo, filename, dictData, data, rules,polarity):
	ColHeadings = ['polarity', 'name','chemical formula','neutral mass', 'rt (min)','protocol', 'file', 'group', 'area', 'height']
	fid = open(filename,'wb')
	for each in ColHeadings:
		fid.write('%s\t' % each)
	fid.write('\n')
	for i,datum in enumerate(data):
		mz = float(dictData[u'compounds'][i][u'mz'])
		mzTol = float(dictData[u'compounds'][i][u'mz_threshold'])
		mzMin = mz - mz*mzTol/1.0e6
		mzMax = mz + mz*mzTol/1.0e6
		rtMin = float(dictData[u'compounds'][i][u'rt_min'])
		rtMax = float(dictData[u'compounds'][i][u'rt_max'])
		cname = dictData[u'compounds'][i]['name']
		cform = dictData[u'compounds'][i]['formula']
		cneum = dictData[u'compounds'][i]['neutral_mass']
		crt = dictData[u'compounds'][i]['rt_peak']
		cprot = "WoM"
		compoundinfo=[cname,cform, cneum, crt, cprot]
		for j,myFile in enumerate(export_fileIds):
			if len(datum)>0:
				tempComp = copy.deepcopy(dictData[u'compounds'][i])
				lowercomp="".join([x.lower() for x in tempComp[u'name']])
				Timea=timeadjustments(tempComp,j,lowercomp,rules,polarity)
				rtMin = float(dictData[u'compounds'][i][u'rt_min'])+Timea
				rtMax = float(dictData[u'compounds'][i][u'rt_max'])+Timea
				idx = np.logical_and( datum[:,2]==myFile, datum[:,0]>=rtMin, datum[:,0]<=rtMax )
				if np.sum(idx)>0:
					x1 = datum[:,0][idx]
					y1 = datum[:,1][idx]
					
					#baseline subtraction of min value between rt_min-extratime and rt_max+extratime
					minY = np.min(datum[:,1][(datum[:,2]==myFile)])
					y2 = y1 - minY  
					#baseline subtraction of min value between rt_min-extratime and rt_max+extratime

					#peakheight
					maxY = np.max(datum[:,1][idx])
					pk_height=maxY-minY
					#peakheight
					
	#                 if myname.startswith('ist'):
	#                     y1 = y1[:]
	#                 else:    
	#                     y1 = y1[:] / fileInfo['normalization_factor'][j]
	#				y2 = y2[:] / fileInfo['normalization_factor'][j]
					parea=np.sum(y2)
				else:
					parea=0
					pk_height=0
			else:
				parea=0
				pk_height=0
			for h,f in enumerate(fileInfo['fid']):
				if f == myFile:
					fpol= fileInfo['polarity'][h]
					fname= fileInfo['name'][h]
					fgroup=fileInfo['group'][h]
			fid.write('%s\t' % fpol)
			for each in compoundinfo:
				fid.write('%s\t' % each)
			fid.write('%s\t%s\t%s\t%s\t\n' % (fname,fgroup,parea,pk_height))
	#     fid.write('\n')
	fid.close()
	print "output file saved to working directory: %s" % filename

def splitresultsforWoM(export_fileIds,resultsfile,fileInfo):
	export_groups = []
	for i,myFile in enumerate(export_fileIds):
		for j,fid in enumerate(fileInfo['fid']):
			if fid == myFile:
				export_groups.append(fileInfo['group'][j])    	
	uniquegroupfiles=[]
	for each in export_groups:
		if '06 - %s.txt' % (each) not in uniquegroupfiles:
			uniquegroupfiles.append('06 - %s.txt' % (each))
	
	with open(resultsfile, 'rU') as results_obj:
		headorder=csv.DictReader(results_obj,dialect='excel-tab')
		ho=headorder.fieldnames    # get the upload order of the fileheadings
		newresults = list(headorder) #dictionary is unordered, use ho to reorder headings and row values

	for each in uniquegroupfiles: 
		print each,
		newfile=open(each,'wb')
		for head in ho[1:]:     #write the first line column names in the original upload order
			newfile.write('%s\t' % head)
		newfile.write('\n')
		for j, line in enumerate(newresults):
			if '%s.txt' %(line['group'])==each:
				newrowtext=[line[i] for i in ho[1:]]
				 # reorders row values into original upload order
				for val in newrowtext:
					newfile.write('%s\t' % val)
				newfile.write('\n')
				#print newrowtext
		newfile.close()
	

#!/usr/bin/env python

'''Implementation of Sliding Window Minimum Algorithm. This module contains
one function `sliding_window_minimum`.
See http://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
for a longer explanation of the algorithm.'''

from collections import deque

def sliding_window_minimum(k, li):
    '''
    A iterator which takes the size of the window, `k`, and an iterable,
    `li`. Then returns an iterator such that the ith element yielded is equal
    to min(list(li)[max(i - k + 1, 0):i+1]).
    Each yield takes amortized O(1) time, and overall the generator takes O(k)
    space.
    __author__ = "Keegan Carruthers-Smith"
	__email__ = "keegan.csmith@gmail.com"
	__license__ = "MIT"
    '''

    window = deque()
    for i, x in enumerate(li):
        while window and window[-1][0] >= x:
            window.pop()
        window.append((x, i))
        while window[0][1] <= i - k:
            window.popleft()
        yield window[0][0]



 
def peakdet(v, delta, x = None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
	
    import sys
    from numpy import NaN, Inf, arange, isscalar, asarray, array
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True
 
    return array(maxtab), array(mintab)

def splitFileListAcrossTwoDimensions(myFiles):
    import re
    fileInfo = {}
    fileInfo['date']=[]
    fileInfo['polarity']=[]
    fileInfo['conc']=[] #this is the concentration of material
    fileInfo['temp']=[] #this is the temperature of the inubation
    fileInfo['group']=[] #this is a unique class of time by group

    for file in myFiles:
        oldName = file
        file = file[:-5]
        file = re.sub('blank\d','blank_0',file)
        file = re.sub('RT','23',file) #substitute the room temperature with 23 for consistency with naming
        fileInfo['polarity'].append(file[-3:])
        file = file[:-4]
        fileInfo['date'].append(file[:6])
        file = file[7:]
        temperature = re.findall(('_[0-9]+$'),file)[0]
        fileInfo['temp'].append(temperature.replace('_',''))
        file = re.sub(temperature+'$','',file)
        file = file[:-2] # strip off the _replicate
        file = file.replace('_','.')
        file = file.replace('bla','blank')
        fileInfo['conc'].append(file)
        fileInfo['group'].append(fileInfo['conc'][-1]+'@'+fileInfo['temp'][-1])
    return fileInfo

def groupFilesAcrossTwoDimensions(fileInfo):
# This block is an example to teach people how to use regular expressions to put files into groups according to their filename
# by naming your files in a consistent manner, it can make analysis of N-way comparisons much quicker
	import re
	fileInfo['date']=[]
	fileInfo['polarity']=[]
	fileInfo['conc']=[] #this is the concentration of material
	fileInfo['temp']=[] #this is the temperature of the inubation
	fileInfo['group']=[] #this is a unique class of time by group

	for file in fileInfo['name']:
	    oldName = file
	    file = file[:-5]
	    file = re.sub('blank\d','blank_0',file)
	    file = re.sub('RT','23',file) #substitute the room temperature with 23 for consistency with naming
	    fileInfo['polarity'].append(file[-3:])
	    file = file[:-4]
	    fileInfo['date'].append(file[:6])
	    file = file[7:]
	    temperature = re.findall(('_[0-9]+$'),file)[0]
	    fileInfo['temp'].append(temperature.replace('_',''))
	    file = re.sub(temperature+'$','',file)
	    file = file[:-2] # strip off the _replicate
	    file = file.replace('_','.')
	    file = file.replace('bla','blank')
	    fileInfo['conc'].append(file)
	    fileInfo['group'].append(fileInfo['conc'][-1]+'@'+fileInfo['temp'][-1])

	fileInfo['pos_groups']={}
	myGroups = np.unique(fileInfo['group'])
	# # print fileInfo['group']
	for group in myGroups:
	    indices = [i for i, elem in enumerate(fileInfo['group']) if group == elem]
	    fileInfo['pos_groups'][group] = []
	    for i in indices:
	        if fileInfo['polarity'][i] == 'pos':
	            fileInfo['pos_groups'][group].append(fileInfo['fid'][i])
	    # print fileInfo['pos_groups'][group]
	# print fileInfo['pos_groups'].keys()

	fileInfo['neg_groups']={}
	myGroups = np.unique(fileInfo['group'])
	# # print fileInfo['group']
	for group in myGroups:
	    indices = [i for i, elem in enumerate(fileInfo['group']) if group == elem]
	    fileInfo['neg_groups'][group] = []
	    for i in indices:
	        if fileInfo['polarity'][i] == 'neg':
	            fileInfo['neg_groups'][group].append(fileInfo['fid'][i])
	    # print group
	    # print fileInfo['neg_groups'][group]
	# print fileInfo['neg_groups'].keys()
	return fileInfo
