from argparse import ArgumentParser
import sys, traceback
import time
import numpy
import getdat
import math
sys.path[:0] = ['/jet/share/lib/python']
import numpy as np
from scipy import interpolate
import scipy
import math
import numpy.fft
import py_flush as pyflush
import readSetParameters
import readGapTreeOld
import getCalibrationFactorSet
import sys
import ppf 
import csv
import load_kk3_raw
import load_ECE
reload(load_kk3_raw)
reload(load_ECE)

#pylab
import numpy
import matplotlib
from matplotlib import pylab, mlab, pyplot
np = numpy
plt = pyplot

from IPython.display import display
from IPython.core.pylabtools import figsize, getfigs

from pylab import *
from numpy import *

def main():
	#~ pulses1 = np.arange(96660,97000)
	#~ pulses4 = np.arange(99360,99999)
	#~ pulses2 = np.arange(97374,98533)
	#~ pulses3 = np.arange(98533,99289)
	#~ pulses = np.concatenate((pulses1,pulses2,pulses3,pulses4))

	pulses_DT_11 = [99815,99817,99818] #hybrid like? Y.Kazakov experiments
	pulses_DT_17 = [99597,99598,99633] #hybrid like
	pulses_DT_09 = [99511,99801,99802,99870,99871] #hybrid like?
	pulses_DT_02 = [99960,99961,99962,99963,99964,99965,99969,99970,99971,99972] #hybrid like, no LIDAR 99530
	pulses_DT_05 = [99981,99982] #hybrid-like
	#pulses_DT_hybridlike = [99815,99817,99818,99597,99598,99633,99801,99802,99981,99982]
	pulses_DT_hybrid = [99448,99449,99450,99452,99455,99541,99542,99543,99544,99594,99595,99596,99760,99761,99866,99867,99868,99869,99908,99910,99912,99914,99949,99950,99951,99953] #no LIDAR 99887,99527
	pulses_DT_baseline = [99512,99520,99795,99796,99797,99799,99805,99861,99862,99863,99878,99943,99944,99948] #99915,99916 look very promising but no LIDAR :(. 99512 has Ne seeding 
	pulses_DT_others = pulses_DT_11+pulses_DT_17+pulses_DT_09+pulses_DT_02+pulses_DT_05

	pulses_DD_baseline_noneon = [96713,96714,96729,96730,96731,96744,96745,96749,96750,96891,96892,96893]
	pulses_DD_hybrid_all = [96951,96953,96954,96955,96957,96958,96765,96767,96768,96769,96770,96771,96773,96776,96777,96779,96717,96718,96719,96720,96721,96722,96723,96724,96725,96726,96663,96664,96665,96668,96670,96671,96674,96676,96677,96678,96679,96680,96681]
	pulses_DD_baseline = [96990,96992,96993,96994,96996,96998,96999]
	pulses_DD_03 = [96838,96840,96848,96850,96851,96852]

	pulses_DD_03_bis = [97771]
	pulses_DD_19_bis = [97511]
	pulses_DD_29_bis = [97382]
	pulses_DD_39_bis = [97484,97486,97490,97493]
	pulses_DD_baseline_bis = [97393,97394,97395,97396,97399,97400,97405,97406,97468,97469,97470,97471,97472,97473,97582,97583,97709,97746,97749,97751,97802,97803,
		97806,97807,97809,97832,97834,97943,97971,97981,97989,97995,98004,98005]
	pulses_DD_hybrid_bis = [97451,97452,97455,97457,97458,97460,97589,97591,97594,97604,97605,97683,97684,97685,97686,97732,97734,97735,97736,97737,
		97740,97741,97742,97779,97780,97781,97783,97784,97785,97786,97787,97788,97790,97791,97796,97797,97798,97799,97844,97847,97852,
		97896,97898,97976,97977]
	pulses_DD_others_bis = pulses_DD_19_bis+pulses_DD_29_bis+pulses_DD_39_bis

	pulses_TT_restart = [98782,98786,98894,98895,98896,98901,98921,98922,98923,98924,98925,99135,99139,99187]	
	pulses_TT_21 = [98791,98885,98891,98892]
	pulses_TT_19 = [98907,98912]
	pulses_TT_hybrid = [98913,99162,99163,99164,99273,99274]
	pulses_TT_10 = [99171]
	pulses_TT_03 = [99194,99195,99198,99206]
	pulses_TT_baseline = [99281,99282]
	pulses_TT_others = pulses_TT_10+pulses_TT_19+pulses_TT_21+pulses_TT_restart

	pulses_all_hybrid = np.array(pulses_DD_hybrid_all+pulses_DD_hybrid_bis+pulses_TT_hybrid+pulses_DT_hybrid)
	pulses_all_baseline = np.array(pulses_DD_baseline+pulses_DD_baseline_noneon+pulses_DD_baseline_bis+pulses_TT_baseline+pulses_DT_baseline)
	pulses_all_03 = np.array(pulses_DD_03+pulses_DD_03_bis+pulses_TT_03)
	pulses_all_others = np.concatenate([pulses_DT_others,pulses_DD_others_bis,pulses_TT_others])
	pulses_all = np.concatenate([pulses_all_hybrid,pulses_all_baseline,pulses_all_03,pulses_all_others])

	pulses = sort(pulses_all)

	pulses = np.array(pulses,dtype=float)
	write_file = 1
	Te_thresh = 1000
	flushid = 'EFIT'
	uidmag = 'jetppf'
	uidK1 = 'mfontana' #'mfontana'
	seqmag = 0
	#~ filename_spec = '_'+str(Te_thresh)+'_'+str(pulses[0])+'-'+str(pulses[-1])+'-'+flushid+'new'
	filename_spec = '_test_1'
	tstep = 0.25 #CHECK IF THIS NEEDS CHANGING		
	readnow = 0

	magyes = 0
	num_samples = 0
	n_k1fail = 0
	n_lidrfail = 0
	n_lowen = 0
	n_lowT = 0

	kB_JperK	= 1.3806504e-23
	e_C			= 1.602176487e-19 	#C
	me_kg		= 9.10938215e-31 	#kg
	mH_kg		= 1.6735575e-27		#kg
	mD_kg		= 3.3435838e-27		#kg
	mT_kg		= 5.0082676e-27		#kg
	epsilon		= 8.854187817e-12 	#F/m
	c0			= 299792458 		#m/s
	mu0			= 1.256637e-6		#N/A^2
	K2eV = kB_JperK/e_C
	c_mpers = 299792458
	constantRayleighJeansLaw_pW = 1e12*kB_JperK/(c_mpers*c_mpers)

	pulseok = []
	pulsenolid = []
	pulsenok1 = []
	pulselowT = []
	pulselowen = []
	pulsenokk1errors = []
	tok = []
	timeidok = []
	lidT = []
	lidTLFS = []
	lidTmax = []
	lidTup = []
	lidTlow = []
	lidne = []
	k1T = []
	k1TLFS = []	
	k1Tmax = []
	ratio = []
	lidT01 = []
	lidT06 = []
	k1T01 = []
	k1T06 = []
	Bt = []
	Ip = []
	Elid = []
	NBIon = []
	ICRHon = []
	Pbeam = []
	PICRH = []
	res1 = []
	res2 = []
	res3 = []
	res4 = []
	res5 = []
	res6 = []
	Tmax2 = []
	Tmax3 = []
	Tmean2 = []
	Tmean3 = []
	Tmean2LFS = []
	Tmean3LFS = []	
	Hconc = []
	Dconc = []
	Tconc = []
	wfast02 = []
	errMean2 = []
	errMean3 = []
	errMax2 = []
	errMax3 = []
	zaxis = []
	betan = []
	gwf = []
	nus = []
	rhos = []
	wdia = []
	q95 = []
	valf = []
	vth = []
	vratio = []

	ratio_filename = open('lid_k1_ratio'+filename_spec+'_zaxis.txt', 'w')

	for ip in range(len(pulses)):
		pulse = pulses[ip]
		print('Processing pulse '+str(pulse))
		Ipla,IplaR,Iplat,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ieripla = ppf.ppfdata(pulse,'MAGN','IPLA',seq=0,reshape=1)
		Bvac,BvacR,Bvact,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierbvac = ppf.ppfdata(pulse,flushid,'BVAC',seq=0,reshape=1)
		zaxistmp,zaxisR,zaxist,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierbvac = ppf.ppfdata(pulse,flushid,'ZMAG',seq=0,reshape=1)	
		if ierbvac!=0:
			#print('No '+str(flushid)+', so using EFIT for pulse '+str(pulse))
			Bvac,BvacR,Bvact,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierbvac2 = ppf.ppfdata(pulse,'EFIT','BVAC',seq=0,reshape=1)			
		## ADIMENSIONAL QUANTITIES
		betanorm,betanR,betant,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierbetan = ppf.ppfdata(pulse,'SCAL','BETN',seq=0,reshape=1)		
		gwfrac,gwfracR,gwfract,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,iergwfrac = ppf.ppfdata(pulse,'SCAL','FGDL',seq=0,reshape=1)	
		nustar,nusR,nust,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,iernustar = ppf.ppfdata(pulse,'SCAL','NUS',seq=0,reshape=1)	
		rhostar,rhosR,rhost,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierrhostar = ppf.ppfdata(pulse,'SCAL','RHOS',seq=0,reshape=1)			
		wpd,wpdR,wpdt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierwpd = ppf.ppfdata(pulse,'SCAL','WPD',seq=0,reshape=1)		
		q95tmp,q95R,q95t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierefit = ppf.ppfdata(pulse,'EFIT','Q95',seq=0,reshape=1)					
		###turn time into number of minutes###
		timeppf,dateppf,ier=ppf.pdstd(pulse)
		if timeppf != 1:
			if len(str(dateppf))==6:
				dayt = float(str(dateppf)[:2])*60*24
				montht = float(str(dateppf)[2:4])*60*24*30
				yeart = float(str(dateppf)[4:6])*60*24*365
			else:
				dayt = float(str(dateppf)[0])*60*24
				montht = float(str(dateppf)[1:3])*60*24*30
				yeart = float(str(dateppf)[3:5])*60*24*365	
			if len(str(timeppf))==6:			
				hourt = float(str(timeppf)[:2])*60
				mint = float(str(timeppf)[2:4])
			else:
				hourt = float(str(timeppf)[0])*60
				mint = float(str(timeppf)[1:3])			
			timeid_temp = mint+hourt+dayt+montht+yeart
			if ip == 0:
				timeid0 = timeid_temp
			timeid = timeid_temp-timeid0
		try: ###LOAD NBI PPF###
			NBIptot,xNBI,tNBI,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierNBI = ppf.ppfdata(pulse,'NBI','PTOT',seq=0,reshape=1) # [W]
			if ierNBI == 0:
				useNBI = 1
			else:
				useNBI = 0	
		except:
			print('Cannot load NBI ppf')	
			useNBI = 0		
		try: ###LOAD ICRH PPF###
			ICRHptot,xICRH,tICRH,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierICRH = ppf.ppfdata(pulse,'ICRH','PTOT',seq=0,reshape=1) # [W]
			if ierICRH == 0:
				useICRH = 1
			else:
				useICRH = 0	
		except:
			print('Cannot load ICRH ppf')	
			useICRH = 0		
		try: ###LOAD H% PPF###
			Hconcppf,xHconc,tHconc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierKS3B = ppf.ppfdata(pulse,'KS3B','HTHD',seq=0,reshape=1) #
			Dconcppf,xHconc,tHconc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierKS3B = ppf.ppfdata(pulse,'KS3B','DTHD',seq=0,reshape=1) #
			Tconcppf,xHconc,tHconc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierKS3B = ppf.ppfdata(pulse,'KS3B','TTTD',seq=0,reshape=1) #
			useKS3B = 0
			useKT5P = 0
			if ierKS3B == 0:
				useKS3B = 1
			else:
				Hconcppf,xHconc,tHconc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierKT5P = ppf.ppfdata(pulse,'KT5P','HTHD',seq=0,reshape=1) # 
				Dconcppf,xHconc,tHconc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierKT5P = ppf.ppfdata(pulse,'KT5P','DTHD',seq=0,reshape=1) # 
				Tconcppf,xHconc,tHconc,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierKT5P = ppf.ppfdata(pulse,'KT5P','TTTD',seq=0,reshape=1) #
				if ierKT5P == 0:
					useKT5P = 1
		except:
			print('Cannot load H conc ppf (KS3B or KT5P)')	
			useKS3B = 0
			useKT5P = 0	
		fastdata=0		
		try: ###LOAD FAST PARTICLES PRESSURE###
			fastdata = 1
			ierpion2 = 1
			ierpencil2 = 1
			wperpion,xwperpion,twperpion,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpion1 = ppf.ppfdata(pulse,'PION','WPER',seq=0,reshape=1,uid='jetppf') # [J/m^-3]
			wparpion,xwparpion,twparpion,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpion1 = ppf.ppfdata(pulse,'PION','WPAR',seq=0,reshape=1,uid='jetppf') # [J/m^-3]
			if ierpion1!=0:
				wperpion,xwperpion,twperpion,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpion2 = ppf.ppfdata(pulse,'PION','WPER',seq=0,reshape=1,uid='chain2')
				wparpion,xwparpion,twparpion,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpion2 = ppf.ppfdata(pulse,'PION','WPAR',seq=0,reshape=1,uid='chain2')
			wperpencil,xwperpencil,twperpencil,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpencil1 = ppf.ppfdata(pulse,'NBP2','WPER',seq=0,reshape=1,uid='jetppf') # [J/m^-3]
			wparpencil,xwparpencil,twparpencil,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpencil1 = ppf.ppfdata(pulse,'NBP2','WPAR',seq=0,reshape=1,uid='jetppf') # [J/m^-3]
			if ierpencil1!=0:
				wperpencil,xwperpencil,twperpencil,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpencil2 = ppf.ppfdata(pulse,'NBP2','WPER',seq=0,reshape=1,uid='chain2')
				wparpencil,xwparpencil,twparpencil,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierpencil2 = ppf.ppfdata(pulse,'NBP2','WPAR',seq=0,reshape=1,uid='chain2')
			if (ierpencil1!=0) & (ierpion1!=0) & (ierpencil2!=0) & (ierpion2!=0):
				print('Cannot load fast particle data (PION or PENCIL)')	
				fastdata=0	
		except:
			print('Cannot load fast particle data (PION or PENCIL) - error')	
			fastdata=0								
		try:
			###LIDAR LOAD DENSITY AND TEMPERATURE###
			r0lidtmp,r0lidR,r0lidt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier1 = ppf.ppfdata(pulse,'LIDR','RMID',seq=0,reshape=1) # [m]
			#print(ier1)
			#print(r0lidtmp)
			telidtmp,telidR,telidt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,iertemp = ppf.ppfdata(pulse,'LIDR','TE',seq=0,reshape=1) # [eV]
			teulidtmp,teulidR,teulidt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','TEU',seq=0,reshape=1) # [eV]
			tellidtmp,tellidR,tellidt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','TEL',seq=0,reshape=1) # [eV]
			nelidtmp,nelidR,nelidt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','NE',seq=0,reshape=1) # [m^-3]
			zlid,zlidR,zlidt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','Z',seq=0,reshape=1) 
			lideng,teng,nwds,titles,units,iereng=getdat.getdat('DC/E3-LA-ENG<MON:004',pulse)
			reslid1tmp,reslid1R,reslid1t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','RES1',seq=0,reshape=1) 
			reslid2tmp,reslid2R,reslid2t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','RES2',seq=0,reshape=1) 
			reslid3tmp,reslid3R,reslid3t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','RES3',seq=0,reshape=1) 
			reslid4tmp,reslid4R,reslid4t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','RES4',seq=0,reshape=1) 
			reslid5tmp,reslid5R,reslid5t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','RES5',seq=0,reshape=1) 
			reslid6tmp,reslid6R,reslid6t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','RES6',seq=0,reshape=1) 
			siglid1tmp,siglid1R,siglid1t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','SIG1',seq=0,reshape=1) 
			siglid2tmp,siglid2R,siglid2t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','SIG2',seq=0,reshape=1) 
			siglid3tmp,siglid3R,siglid3t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','SIG3',seq=0,reshape=1) 
			siglid4tmp,siglid4R,siglid4t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','SIG4',seq=0,reshape=1) 
			siglid5tmp,siglid5R,siglid5t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','SIG5',seq=0,reshape=1) 
			siglid6tmp,siglid6R,siglid6t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'LIDR','SIG6',seq=0,reshape=1) 	
			if (mean(r0lidtmp)==0) or (nt==1):
				iertemp = 100	
			if iertemp!=0:
				print('No LIDAR data for pulse '+str(pulse))
				n_lidrfail = n_lidrfail+1
				pulsenolid = pulsenolid+[pulse]		
				continue								
			dlideng = diff(lideng)
			tlideng_jump = teng[dlideng>0.05]
			r0lidtmp[r0lidtmp==0] = nan
			telidtmp[telidtmp==0] = nan
			teulidtmp[teulidtmp==0] = nan
			tellidtmp[tellidtmp==0] = nan
			nelidtmp[telidtmp==0] = nan			
		
		except:
			print('Issues with loading LIDR ppf data for pulse '+str(pulse))
			n_lidrfail = n_lidrfail+1
			pulsenolid = pulsenolid+[pulse]	
			continue
		try:
			k1T1,k1R1,k1t,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'ECM1','PRFL',reshape=1,uid=uidK1,seq=0)
			#k1T,k1R,k1spt,fluxint = load_ECE.load_kk1_equi_simple(pulse,uidK1,flushid,uidmag,seqmag)
			k1sp,k1spfreq,k1spt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ier = ppf.ppfdata(pulse,'ECM1','CSPC',reshape=1,uid=uidK1,seq=0)	
			fmax,fmaxfreq,fmaxt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierkk1 = ppf.ppfdata(pulse,'ECM1','FMAX',reshape=1,uid=uidK1)	
			fcom,fcomfreq,fcomt,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,ierkk1 = ppf.ppfdata(pulse,'ECM1','FCOM',reshape=1,uid=uidK1)	
			try:
				cser,cserfreq,csert,nd,nx,nt,dunits,xunits,tunits,desc,comm,seq,iercser = ppf.ppfdata(pulse,'ECM1','CSER',reshape=1,uid=uidK1)	
				if iercser!=0:
					print('No errors saved for pulse '+str(pulse))
					pulsenokk1errors = pulsenokk1errors+[pulse]
					cser = np.zeros(shape(k1sp))
				cserrel = cser/k1sp	
			except:
				print('Error - No errors saved for pulse '+str(pulse))
				pulsenokk1errors = pulsenokk1errors+[pulse]
				cser = np.zeros(shape(k1sp))
				cserrel = cser/k1sp	
		except:
			print('Issues with KK1 data for pulse '+str(pulse))	
			n_k1fail = n_k1fail+1
			pulsenok1 = pulsenok1+[pulse]
			continue	
		if ierkk1!=0:
			print('Issues with KK1 data for pulse '+str(pulse))	
			pulsenok1 = pulsenok1+[pulse]
			continue				
		###Cycle over the LIDR times, look at the core values###
		tnow = r0lidt[0]
		while tnow<r0lidt[-1]+0.2:
			lidT01now = nan
			lidT06now = nan
			lidTupnow = nan
			lidTlownow = nan
			lidTnow = nan
			lidTnowLFS = nan
			lidnenow = nan
			res1now = nan
			res2now = nan
			res3now = nan
			res4now = nan
			res5now = nan	
			res6now = nan																																						
			engnow = nan
			if min(abs(tnow-r0lidt))<0.2:
				it = np.argmin(abs(tnow-r0lidt))
				tnow = r0lidt[it]
				iteng_jump = np.argmin(abs(tlideng_jump-tnow))
				teng_jump = tlideng_jump[iteng_jump]
				iteng = np.argmin(abs(teng-teng_jump))
				engnow = mean(lideng[iteng+1:iteng+6])
				r0lid = r0lidtmp[it,:]
				telid = telidtmp[it,:]
				idmaxLID = np.nanargmax(abs(telid))
				if isnan(idmaxLID)==0:
					idmaxLID1 = np.nanargmin(abs(r0lid-(r0lid[idmaxLID]-0.1)))
					idmaxLID2 = np.nanargmin(abs(r0lid-(r0lid[idmaxLID]+0.1)))
					lidTmaxnow = np.mean(telid[idmaxLID1:idmaxLID2][np.isnan(telid[idmaxLID1:idmaxLID2])==0])	
				else:
					lidTmaxnow = nan					
				id1 = np.nanargmin(abs(r0lid-2.85))
				id2 = np.nanargmin(abs(r0lid-3.15))
				id3 = np.nanargmin(abs(r0lid-3.0))
				id4 = np.nanargmin(abs(r0lid-3.15))
				lidTnow = np.mean(telid[id1:id2][np.isnan(telid[id1:id2])==0])	
				lidTnowLFS = np.mean(telid[id3:id4][np.isnan(telid[id3:id4])==0])	
				if (engnow>0.4) or (lidTnow>Te_thresh):
					teulid = teulidtmp[it,:]
					tellid = tellidtmp[it,:]
					nelid = nelidtmp[it,:]
					reslid1 = reslid1tmp[it,:]
					reslid2 = reslid2tmp[it,:]
					reslid3 = reslid3tmp[it,:]
					reslid4 = reslid4tmp[it,:]
					reslid5 = reslid5tmp[it,:]
					reslid6 = reslid6tmp[it,:]
					siglid1 = siglid1tmp[it,:]
					siglid2 = siglid2tmp[it,:]
					siglid3 = siglid3tmp[it,:]
					siglid4 = siglid4tmp[it,:]
					siglid5 = siglid5tmp[it,:]
					siglid6 = siglid6tmp[it,:]				
					lidnenow = np.mean(nelid[id1:id2][np.isnan(nelid[id1:id2])==0])/1e19
					res1now = np.mean(reslid1[id1:id2][np.isnan(telid[id1:id2])==0])/np.mean(siglid1[id1:id2][np.isnan(telid[id1:id2])==0])
					res2now = np.mean(reslid2[id1:id2][np.isnan(telid[id1:id2])==0])/np.mean(siglid2[id1:id2][np.isnan(telid[id1:id2])==0])
					res3now = np.mean(reslid3[id1:id2][np.isnan(telid[id1:id2])==0])/np.mean(siglid3[id1:id2][np.isnan(telid[id1:id2])==0])
					res4now = np.mean(reslid4[id1:id2][np.isnan(telid[id1:id2])==0])/np.mean(siglid4[id1:id2][np.isnan(telid[id1:id2])==0])
					res5now = np.mean(reslid5[id1:id2][np.isnan(telid[id1:id2])==0])/np.mean(siglid5[id1:id2][np.isnan(telid[id1:id2])==0])
					res6now = np.mean(reslid6[id1:id2][np.isnan(telid[id1:id2])==0])/np.mean(siglid6[id1:id2][np.isnan(telid[id1:id2])==0])	
					lidTupnow = abs(np.mean(teulid[id1:id2][np.isnan(teulid[id1:id2])==0])-lidTnow)/lidTnow
					lidTlownow = abs(np.mean(tellid[id1:id2][np.isnan(tellid[id1:id2])==0])-lidTnow)/lidTnow
					#lidTupnow = np.sqrt(sum(teulid[id1:id2][np.isnan(teulid[id1:id2])==0]**2))/(id2-id1-1)
					#lidTlownow = np.sqrt(sum(tellid[id1:id2][np.isnan(tellid[id1:id2])==0]**2))/(id2-id1-1)		
					if magyes == 1:
						pyflush.flushinit(15, int(pulse), r0lidt[it], 0, seqmag, uidmag, flushid, 0)
						spatialLimits = []
						spatialLimits = pyflush.Flush_getSplineDomainLimits()
						flux = pyflush.Flush_getFlux(len(telidR),numpy.float64(telidR)*100,numpy.float64(zlid)*100)		
						lidflux = numpy.sqrt(flux[0])
						idLFSlid = np.nanargmin(lidflux)
						telidLFS = telid[idLFSlid:]
						lidTnow01 = np.mean(telid[(lidflux<0.12)&(np.isnan(telid)==0)])		
						lidTnow06 = np.mean(telidLFS[(lidflux[idLFSlid:]<0.4)&(lidflux[idLFSlid:]>0.3)&(np.isnan(telidLFS)==0)])	
						#lidTnow06 = np.mean(telidLFS[(lidflux[idLFSlid:]<0.67)&(lidflux[idLFSlid:]>0.53)&(np.isnan(telidLFS)==0)])	
						lidT01 = lidT01+[lidTnow01]
						lidT06 = lidT06+[lidTnow06]
				#~ elif lidTnow<Te_thresh:
						#~ #print('Te LIDAR = '+str(lidTnow)+' for pulse '+str(pulse))	
						#~ n_lowT = n_lowT+1
						#~ pulselowT = pulselowT+[pulse]
				#~ elif engnow<0.4:
						#~ #print('Laser energy = '+str(engnow)+' for pulse '+str(pulse))	
						#~ n_lowen = n_lowen+1	
						#~ pulselowen = pulselowen+[pulse]		
								
			#Process KK1 and other quantities
			itbvac = np.argmin(abs(Bvact-tnow))
			Btnow = Bvac[itbvac]	
			itipla = np.argmin(abs(Iplat-tnow))
			Ipnow = Ipla[itipla]	
			itzaxis = np.argmin(abs(zaxist-tnow))
			zaxisnow = zaxistmp[itzaxis]
			itq95 = np.argmin(abs(q95t-tnow))
			q95now = q95tmp[itq95]
			if ierbetan == 0:
				itscal = np.argmin(abs(betant-tnow))
				betanormnow = betanorm[itscal]
				gwfracnow = gwfrac[itscal]
				nustarnow = nustar[itscal]
				rhostarnow = rhostar[itscal]		
				wpdnow = wpd[itscal]
			else:
				betanormnow = nan
				gwfracnow = nan
				nustarnow = nan
				rhostarnow = nan	
				wpdnow = nan					
			itk1 = numpy.argmin(abs(k1t-tnow))			
			if abs(k1t[itk1]-tnow)<0.015:
				radiationTemperature = np.zeros((len(k1spfreq)))
				for i in range(len(k1spfreq)):
					radiationTemperature[i] = numpy.multiply(k1sp[itk1,i],1/(constantRayleighJeansLaw_pW*k1spfreq[i]**2))*K2eV #[eV]
				Rs = numpy.float64(k1R1)
				zkk1 = 0.356 #m	
				zkk1LOS = np.ones([len(Rs)])*zkk1
				idmaxECE = np.nanargmax(abs(k1T1[itk1,:]))
				if isnan(idmaxECE)==0:
					idmaxECE1 = np.nanargmin(abs(k1R1-(k1R1[idmaxECE]-0.1)))				
					idmaxECE2 = np.nanargmin(abs(k1R1-(k1R1[idmaxECE]+0.1)))				
					k1Tmaxnow = np.mean(k1T1[itk1,idmaxECE1:idmaxECE2])	
				else:
					k1Tmaxnow = nan					
				idk11 = np.argmin(abs(k1R1-2.85))
				idk12 = np.argmin(abs(k1R1-3.15))		
				idk13 = np.argmin(abs(k1R1-3.0))
				idk14 = np.argmin(abs(k1R1-3.15))	
				k1Tnow = np.mean(k1T1[itk1,idk11:idk12])
				k1TLFSnow = np.mean(k1T1[itk1,idk13:idk14])				
				idfmax2 = np.argmin(abs(fmax[itk1]*2-k1spfreq))
				idfmax3 = np.argmin(abs(fmax[itk1]*2.9-k1spfreq)) #2.9 instead of 3 because of G.Giruzzi comment
				#~ Tmax2now = radiationTemperature[idfmax2]
				#~ Tmax3now = radiationTemperature[idfmax3]
				idfcom2 = np.argmin(abs(fcom[itk1]*2-k1spfreq))
				idfcom3 = np.argmin(abs(fcom[itk1]*2.9-k1spfreq))				
				idf5minus2 = np.argmin(abs(fcom[itk1]*2-fcom[itk1]*2*0.05-k1spfreq))
				idf5plus2 = np.argmin(abs(fcom[itk1]*2+fcom[itk1]*2*0.05-k1spfreq))+1
				idf5minus3 = np.argmin(abs(fcom[itk1]*2.9-fcom[itk1]*2.9*0.05-k1spfreq))
				idf5plus3 = np.argmin(abs(fcom[itk1]*2.9+fcom[itk1]*2.9*0.05-k1spfreq))+1
				errMean2now = np.sqrt(sum(cserrel[itk1,idf5minus2:idf5plus2]**2))/(idf5plus2-idf5minus2-1)
				errMean3now = np.sqrt(sum(cserrel[itk1,idf5minus3:idf5plus3]**2))/(idf5plus3-idf5minus3-1)									
				Tmean2now = np.mean(radiationTemperature[idf5minus2:idf5plus2])
				Tmean3now = np.mean(radiationTemperature[idf5minus3:idf5plus3])			
				Tmean2LFSnow = np.mean(radiationTemperature[idf5minus2:idfcom2])
				Tmean3LFSnow = np.mean(radiationTemperature[idf5minus3:idfcom3])						
				if (radiationTemperature[idf5minus2]>radiationTemperature[idfcom2])&(radiationTemperature[idf5plus2]>radiationTemperature[idfcom2]):
					Tmax2now = radiationTemperature[idfcom2]
					Tmax3now = radiationTemperature[idfcom3]
					errMax2now = cserrel[itk1,idfcom2]
					errMax3now = cserrel[itk1,idfcom3]					
				else:
					Tmax2now = np.max(radiationTemperature[idf5minus2:idf5plus2])
					Tmax3now = np.max(radiationTemperature[idf5minus3:idf5plus3])
					idmax2 = np.argmax(radiationTemperature[idf5minus2:idf5plus2])
					idmax3 = np.argmax(radiationTemperature[idf5minus3:idf5plus3])
					errMax2now = cserrel[itk1,idf5minus2:idf5plus2][idmax2]
					errMax3now = cserrel[itk1,idf5minus3:idf5plus3][idmax3]					
				if magyes == 1:
					pyflush.flushinit(15, int(pulse), k1t[itk1], 0, seqmag, uidmag, flushid, 0)
					BTor = pyflush.Flush_getBt(len(Rs),Rs*100,zkk1LOS*100)
					Bz   = pyflush.Flush_getBz(len(Rs),Rs*100,zkk1LOS*100)
					Br   = pyflush.Flush_getBr(len(Rs),Rs*100,zkk1LOS*100)
					BLOS = numpy.zeros(len(Rs))
					BLOS = (BTor[0][:]**2+Bz[0][:]**2+Br[0][:]**2)**0.5/1e4 #[T]
					flux = pyflush.Flush_getFlux(len(Rs),Rs*100,zkk1LOS*100)		
					k1flux = numpy.sqrt(flux[0])
					freqs = 2*e_C*BLOS/me_kg/2/3.1415
					k1Tint = np.interp(freqs/1e9,k1spfreq/1e9,radiationTemperature)							
					idLFSk1 = np.nanargmin(k1flux)
					tek1LFS = k1Tint[idLFSk1:]
					k1Tnow01 = np.mean(k1Tint[(k1flux<0.12)&(np.isnan(k1Tint)==0)])		
					k1Tnow06 = np.mean(tek1LFS[(k1flux[idLFSk1:]<0.4)&(k1flux[idLFSk1:]>0.3)&(np.isnan(tek1LFS)==0)])	
					#k1Tnow06 = np.mean(tek1LFS[(k1flux[idLFSk1:]<0.67)&(k1flux[idLFSk1:]>0.53)&(np.isnan(tek1LFS)==0)])		
				if k1Tnow!=0:
					pulseok = pulseok+[pulse]
					tok = tok+[r0lidt[it]]
					timeidok = timeidok+[timeid]
					k1Tmax = k1Tmax+[k1Tmaxnow]
					k1T = k1T+[k1Tnow]
					k1TLFS = k1TLFS+[k1TLFSnow]
					Tmax2 = Tmax2+[Tmax2now]
					Tmax3 = Tmax3+[Tmax3now]
					Tmean2 = Tmean2+[Tmean2now]
					Tmean3 = Tmean3+[Tmean3now]		
					Tmean2LFS = Tmean2LFS+[Tmean2LFSnow]
					Tmean3LFS = Tmean3LFS+[Tmean3LFSnow]
					errMean2 = errMean2+[errMean2now]
					errMean3 = errMean3+[errMean3now]	
					errMax2 = errMax2+[errMax2now]
					errMax3 = errMax3+[errMax3now]								
					Bt = Bt+[Btnow]	
					Ip = Ip+[Ipnow]		
					zaxis = zaxis+[zaxisnow]
					lidTmax = lidTmax+[lidTmaxnow]		
					lidT = lidT+[lidTnow]
					lidTLFS = lidTLFS+[lidTnowLFS]
					lidTup = lidTup+[lidTupnow]
					lidTlow = lidTlow+[lidTlownow]
					lidne = lidne+[lidnenow]
					res1 = res1+[res1now]
					res2 = res2+[res2now]
					res3 = res3+[res3now]
					res4 = res4+[res4now]
					res5 = res5+[res5now]		
					res6 = res6+[res6now]																																						
					Elid = Elid+[engnow]		
					betan = betan+[betanormnow]	
					gwf = gwf+[gwfracnow]
					nus = nus+[nustarnow]
					rhos = rhos+[rhostarnow]
					wdia = wdia+[wpdnow]
					q95 = q95+[q95now]									
					num_samples = num_samples+1
					if magyes == 1:
						lidT01 = lidT01+[lidTnow01]
						lidT06 = lidT06+[lidTnow06]
						k1T01 = k1T01+[k1Tnow01]
						k1T06 = k1T06+[k1Tnow06]
					if useICRH == 1:
						iticrh = np.argmin(abs(tICRH-tnow))
						if ICRHptot[iticrh]>0.5e6:
							ICRHon = ICRHon+[1]
							PICRHnow = ICRHptot[iticrh]/1e6
						else:
							ICRHon = ICRHon+[0]	
							PICRHnow = 0
					else:
						ICRHon = ICRHon+[0]	
						PICRHnow = 0	
					PICRH = PICRH+[PICRHnow]						
					if useNBI == 1:
						itnbi = np.argmin(abs(tNBI-tnow))
						if NBIptot[itnbi]>1e6:
							NBIon = NBIon+[1]
							Pbeamnow = NBIptot[itnbi]/1e6
						else:
							NBIon = NBIon+[0]	
							Pbeamnow = 0	
					else:
						NBIon = NBIon+[0]	
						Pbeamnow = 0	
					Pbeam = Pbeam+[Pbeamnow]
					if useKT5P==1 or useKS3B==1:
						ithconc = np.argmin(abs(tHconc-tnow))		
						Hconcnow = Hconcppf[ithconc]	
						Dconcnow = Dconcppf[ithconc]	
						Tconcnow = Tconcppf[ithconc]	
						valfnow = -Btnow/sqrt(mu0*lidnenow*1e19*(Hconcnow*mH_kg+Dconcnow*mD_kg+Tconcnow*mT_kg))
						vthnow = sqrt(lidTnow*e_C/((Hconcnow*mH_kg+Dconcnow*mD_kg+Tconcnow*mT_kg)))
					else:
						Hconcnow = 0
						Dconcnow = 0
						Tconcnow = 0	
						valfnow = -Btnow/sqrt(mu0*lidnenow*1e19*mD_kg)
						vthnow = sqrt(lidTnow*e_C/mD_kg)						
					Hconc = Hconc+[Hconcnow]
					Dconc = Dconc+[Dconcnow]
					Tconc = Tconc+[Tconcnow]
					valf = valf+[valfnow]
					vth = vth+[vthnow]
					vratio = vratio+[valfnow/vthnow]
					if fastdata==1:
						wfastmean = nan
						if ICRHon[-1]==1 or ((ierpencil1!=0) and (ierpencil2!=0)):
							if (ierpion1==0) or (ierpion2==0):
								itwper = np.argmin(abs(twperpion-tnow))
								iwper02 = np.argmin(abs(xwperpion-0.2))
								if itwper != 0:
									wfastmean = mean(wperpion[itwper,:iwper02+1]+wparpion[itwper,:iwper02+1])
								else:
									wfastmean = mean(wperpion[:iwper02+1]+wparpion[:iwper02+1])	
								
						else:
							if (ierpencil1==0) or (ierpencil2==0):
								itwper = np.argmin(abs(twperpencil-tnow))
								iwper02 = np.argmin(abs(xwperpencil-0.2))
								if itwper != 0:
									wfastmean = mean(wperpencil[itwper,:iwper02+1]+wparpencil[itwper,:iwper02+1])
								else:	
									wfastmean = mean(wperpencil[:iwper02+1]+wparpencil[:iwper02+1])		
						if (~isnan(wfastmean))&(Pbeamnow==0)&(PICRHnow==0):
							wfastmean = nan					
					else:
						wfastmean = nan
					wfast02 = wfast02+[wfastmean]	
					###### WRITE FILES FOR PROFILES ######
					if write_file==1:
						print >> ratio_filename, '%6d %4.3f %8.2f %5.3f %4.3f %4.2f %2d %2d %6.6f %6.6f\
						 %6.6f %6.6f %6.6f %6.6f %4.3f %4.2f %5.2f %6.2f %6.2f %6.2f\
						 %6.2f %6.6f %6.6f %6.6f %6.6f %4.2f %4.2f %6.2f %5.3f %5.3f \
						 %8.2f %4.2f %8.2f %4.2f %5.3f %8.2f %8.2f\
						 %5.3f %5.3f %8.6f %8.6f %5.3f %8.2f %8.2f %5.3f'\
						  %(pulse,r0lidt[it],lidTnow,lidTnow/k1Tnow,engnow,Btnow,NBIon[-1],ICRHon[-1],res1now,res2now,\
						  res3now,res4now,res5now,res6now,lidnenow,Ipnow/1e6,Pbeamnow,Tmax2now,Tmax3now,Tmean2now,\
						  Tmean3now,errMean2now,errMean3now,errMax2now,errMax3now,PICRHnow,Hconcnow*100,wfastmean,lidTupnow,lidTlownow,\
						  lidTnowLFS,lidTnowLFS/k1TLFSnow,lidTmaxnow,lidTmaxnow/k1Tmaxnow,zaxisnow,Tmean2LFSnow,Tmean3LFSnow,\
						  betanormnow,gwfracnow,nustarnow,rhostarnow,q95now,valfnow,vthnow,valfnow/vthnow)																						
			tnow = tnow + tstep					

	ratio_filename.close()

try:
    main()
except SystemExit:
    traceback.print_exc(file=sys.stderr)
    raise
except:
    print "-------unhandled exception!---------"
    traceback.print_exc(file=sys.stderr)
    sys.exit(127)
