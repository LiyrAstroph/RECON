#/*
# * RECON Copyright (C) 2018 Yan-Rong Li
# * A package for measuring spectral power and reconstructing time series in AGN.
# * 
# * Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# * 
# * implement Bayesian posterior predictive checking.
# */

import os
import copy
import numpy as np
import numpy.fft as fft
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming import pyPDM
from astropy.stats import LombScargle

#=======================================================
# resmple time series
#
#=======================================================
def resample(t, f):
  trs = np.linspace(t[0], t[-1], len(t))
  frs = np.interp(trs, t, f)
  return trs, frs

#=======================================================
# calculate PSD
#
#=======================================================
def fft_psd(t, f, normed=True):
  
  # interpolate light cruve into evenly sampling
  trs, frs = resample(t, f)
  
  mean = np.mean(frs)
  
  # set the starting time to be zero and remove the mean
  #ts = ts - ts[0]
  frs = frs- mean
  
  # fft
  pfft = fft.rfft(frs) 

  # frequency, the positive part
  freq = np.arange(1, len(trs)/2+1).astype(float)/(trs[-1]-trs[0])

  # power spectral density
  psd = np.absolute(np.real(pfft[1:]))**2 + np.absolute(np.imag(pfft[1:]))**2
  
  # normalization
  if normed == True:
    psd *= 2.0 * (trs[-1]-trs[0])/(len(trs)**2)

  return freq, psd

#=============================================================
# rebin PSD
#
#=============================================================
def fft_psd_rebin(freq, psd):
  freq_work = np.log10(freq)
  psd_work = np.log10(psd)

  freq0 = freq[0]
  width = freq[1] - freq[0]
  width_log = np.log10(1.1)
  nbins_min = 1

  ic = 0;
  f0 = freq_work[0];
  i0 = 0

  freq_rebin = np.zeros(len(freq))
  psd_rebin = np.zeros(len(freq))

  while f0 < freq_work[-1]:
    f1 = f0 + width_log

    i1 = int( (10.0**f1 - freq0)/width )

    if i1-i0 < nbins_min:
      i1 = i0 + nbins_min

    if ((len(freq)-1 - i1 < nbins_min)) or (freq_work[-1] - f1 < width_log):
      i1 = len(freq) - 1

    f1 = freq_work[i1]
    freq_rebin[ic] = 0.5*(freq_work[i0] + freq_work[i1])
    fsum = freq_work[i1] - freq_work[i0];
    psum = 0.0
    for i in range(i0, i1):
      psum += 0.5*(psd_work[i] + psd_work[i+1]) * (freq_work[i+1] - freq_work[i])

    psd_rebin[ic] = psum/fsum

    i0 = i1
    f0 = f1
    ic += 1

  return 10.0**freq_rebin[:ic], 10.0**psd_rebin[:ic]

#=======================================================
# single power law PSD function
#
#=======================================================
def psd_power_law(fk, arg):
  A= np.exp(arg[0])
  alpha=arg[1]
  cnoise = np.exp(arg[2])
  
  psd = np.zeros(fk.shape)
  
  idx = np.where(fk >= freq_limit)
  psd[idx[0]] = A * fk[idx[0]]**(-alpha) + cnoise
  idx = np.where(fk < freq_limit)
  psd[idx[0]] = A * freq_limit**(-alpha) + cnoise
  
  return psd

#=======================================================
# square root of single power law PSD function
#
#=======================================================
def psd_power_law_sqrt(fk, arg):
  A= np.exp(arg[0])
  alpha=arg[1]
  cnoise = np.exp(arg[2])
  
  psd = np.zeros(fk.shape)
  
  idx = np.where(fk >= freq_limit)
  psd[idx[0]] = A * fk[idx[0]]**(-alpha) + cnoise
  idx = np.where(fk < freq_limit)
  psd[idx[0]] = A * freq_limit**(-alpha) + cnoise
  
  return np.sqrt(psd)

#=======================================================
# Lorentzian PSD function
#
#=======================================================
def psd_period_lorentz(fk, arg):
  Ap = np.exp(arg[0])
  center = np.exp(arg[1])
  width = np.exp(arg[2])
  psd = Ap/np.pi * width/(width*width + (fk-center)**2)
  
  return psd

#=======================================================
# square root of Lorentzian PSD function  
#
#=======================================================
def psd_period_sqrt_lorentz(fk, arg):
  Ap = np.exp(arg[0])
  center = np.exp(arg[1])
  width = np.exp(arg[2])
  psd = Ap/np.pi * width/(width*width + (fk-center)**2)
  
  return np.sqrt(psd)

#=======================================================
# Gaussian PSD function
#
#=======================================================
def psd_period_gaussian(fk, arg):
  Ap = np.exp(arg[0])
  center = np.exp(arg[1])
  width = np.exp(arg[2])
  psd = Ap/np.sqrt(2.0*np.pi)/width * np.exp( -0.5*(fk-center)**2/width**2 )
  
  return psd

#=======================================================
# square root of Gaussian PSD function  
#
#=======================================================
def psd_period_sqrt_gaussian(fk, arg):
  Ap = np.exp(arg[0])
  center = np.exp(arg[1])
  width = np.exp(arg[2])
  psd = Ap/np.sqrt(2.0*np.pi)/width * np.exp( -0.5*(fk-center)**2/width**2 )
  
  return np.sqrt(psd)

#=======================================================
# generate time series  
#
#=======================================================
def genlc(model):
  global num_params_psd, num_params_psd_sto, num_params_psd_per, nd_sim
  global DT, flux_scale, flux_mean
  global psd_period, psd_period_sqrt
  
  arg = model[:num_params_psd_sto]
  
  fft_work = np.zeros(nd_sim/2+1, dtype=complex)
  
  fft_work[0] = model[num_params_psd]+1j*0.0
  
  freq = 1.0/(nd_sim * DT) * np.linspace(0.0, nd_sim/2, nd_sim/2+1)
  fft_work[1:nd_sim/2] = psd_power_law_sqrt(freq[1:nd_sim/2], arg)/np.sqrt(2.0) \
                        * (model[num_params_psd + 1:num_params_psd + nd_sim -1:2] + 1j*model[num_params_psd + 2:num_params_psd +nd_sim-1:2])
  fft_work[nd_sim/2] = psd_power_law_sqrt(freq[nd_sim/2:], arg) * (model[num_params_psd + nd_sim -1] + 1j*0.0)
  
  if num_params_psd_per > 0:
    arg = model[num_params_psd_sto:num_params_psd]
    fft_work[1:] += psd_period_sqrt(freq[1:], arg) * np.exp(1j*model[num_params_psd + nd_sim:]*2.0*np.pi)
  
  fs = fft.irfft(fft_work) * nd_sim # note the factor 1/n in numpy ifft,
  
  norm = 1.0/np.sqrt(nd_sim) * np.sqrt(nd_sim/(2.0*nd_sim * DT))
  
  ts = DT * ( np.linspace(0, nd_sim-1, nd_sim) - nd_sim/2.0) + time_media
  fs = fs*norm*flux_scale + flux_mean
  
  return ts, fs 

#=======================================================
# generate time series on the time grid of data  
# note that Gaussian noises for the measurement error are included.
#
#=======================================================
def genlc_data(model):
  global num_params_psd, num_params_psd_sto, num_params_psd_per, nd_sim
  global DT, flux_scale, flux_mean
  global lc
  global psd_period, psd_period_sqrt
  
  arg = model[:num_params_psd_sto]
  
  fft_work = np.zeros(nd_sim/2+1, dtype=complex)
  
  fft_work[0] = model[num_params_psd]+1j*0.0
  
  freq = 1.0/(nd_sim * DT) * np.linspace(0.0, nd_sim/2, nd_sim/2+1)
  fft_work[1:nd_sim/2] = psd_power_law_sqrt(freq[1:nd_sim/2], arg)/np.sqrt(2.0) \
                        * (model[num_params_psd + 1:num_params_psd + nd_sim -1:2] + 1j*model[num_params_psd + 2:num_params_psd +nd_sim-1:2])
  fft_work[nd_sim/2] = psd_power_law_sqrt(freq[nd_sim/2:], arg) * (model[num_params_psd + nd_sim - 1]  + 1j*0.0)
  if num_params_psd_per > 0:
    arg = model[num_params_psd_sto:num_params_psd]
    fft_work[1:] += psd_period_sqrt(freq[1:], arg) * np.exp(1j*model[num_params_psd + nd_sim:]*2.0*np.pi)
  
  fs = fft.irfft(fft_work) * nd_sim # note the factor 1/n in numpy ifft,
  
  norm = 1.0/np.sqrt(nd_sim) * np.sqrt(nd_sim/(2.0*nd_sim * DT))
  
  ts = DT * ( np.linspace(0, nd_sim-1, nd_sim) - nd_sim/2.0) + time_media
  fs = fs*norm*flux_scale + flux_mean
  
  fsd = np.interp(lc[:, 0], ts, fs)
  
  fsd += np.random.randn(lc.shape[0]) * ls[:, 2]
  return lc[:, 0], fsd 

#=======================================================  
# generate time series with given PSD
#
#=======================================================
def genlc_psd(model):
  global num_params_psd, num_params_psd_sto, num_params_psd_per
  global DT, flux_scale, flux_mean
  global psd_period, psd_period_sqrt
  
  nd_sim = 1000
  W = 4
  
  arg = model[:num_params_psd_sto]
  
  fft_work = np.zeros(nd_sim/2+1, dtype=complex)
  
  fft_work[0] = np.random.randn()+1j*0.0
  
  freq = 1.0/(nd_sim * DT) * np.linspace(0.0, nd_sim/2, nd_sim/2+1)
  fft_work[1:nd_sim/2] = psd_power_law_sqrt(freq[1:nd_sim/2], arg)/np.sqrt(2.0) \
                        * (np.random.randn(nd_sim/2-1) + 1j*np.random.randn(nd_sim/2-1))
  fft_work[nd_sim/2] = psd_power_law_sqrt(freq[nd_sim/2:], arg) * (np.random.randn() + 1j*0.0)
  
  if num_params_psd_per > 0:
    arg = model[num_params_psd_sto:num_params_psd]
    fft_work[1:] += psd_period_sqrt(freq[1:], arg) * np.exp(1j*np.random.rand(nd_sim/2)*2.0*np.pi)
  
  fs = fft.irfft(fft_work) * nd_sim # note the factor 1/n in numpy ifft,
  
  norm = 1.0/np.sqrt(nd_sim) * np.sqrt(nd_sim/(2.0*nd_sim * DT))
  
  ts = DT * ( np.linspace(0, nd_sim-1, nd_sim) - nd_sim/2.0) + time_media
  fs = fs*norm*flux_scale + flux_mean
  
  return ts, fs 

#=======================================================
# generate time series on the observed time grid with given PSD
# note that Gaussian noises for the measurement error are included.
#
#=======================================================
def genlc_psd_data(model):
  global num_params_psd, num_params_psd_sto, num_params_psd_per, nd_sim
  global DT, flux_scale, flux_mean
  global lc
  global psd_period, psd_period_sqrt
  
  arg = model[:num_params_psd_sto]
  
  fft_work = np.zeros(nd_sim/2+1, dtype=complex)
  
  fft_work[0] = np.random.randn()+1j*0.0
  
  freq = 1.0/(nd_sim * DT) * np.linspace(0.0, nd_sim/2, nd_sim/2+1)
  fft_work[1:nd_sim/2] = psd_power_law_sqrt(freq[1:nd_sim/2], arg)/np.sqrt(2.0) \
                        * (np.random.randn(nd_sim/2-1) + 1j*np.random.randn(nd_sim/2-1))
  fft_work[nd_sim/2] = psd_power_law_sqrt(freq[nd_sim/2:], arg) * (np.random.randn() + 1j*0.0)
  
  if num_params_psd_per > 0:
    arg = model[num_params_psd_sto:num_params_psd]
    fft_work[1:] += psd_period_sqrt(freq[1:], arg) * np.exp(1j*np.random.rand(nd_sim/2)*2.0*np.pi)
  
  fs = fft.irfft(fft_work) * nd_sim # note the factor 1/n in numpy ifft,
  
  norm = 1.0/np.sqrt(nd_sim) * np.sqrt(nd_sim/(2.0*nd_sim * DT))
  
  ts = DT * ( np.linspace(0, nd_sim-1, nd_sim) - nd_sim/2.0) + time_media
  fs = fs*norm*flux_scale + flux_mean
  
  fsd = np.interp(lc[:, 0], ts, fs)
  
  fsd += np.random.randn(lc.shape[0]) * lc[:, 2]
  return lc[:, 0], fsd 

#=======================================================
# load posterior sample
#
#=======================================================
def load_sample():
  global sample
  sample = np.loadtxt("data/posterior_sample.txt")

  print "AR:", np.mean(sample[:, 0]/np.log(10.0)) + 2.0*np.log10(flux_scale), np.std(sample[:, 0]/np.log(10.0))
  print "alpha:", np.mean(sample[:, 1]), np.std(sample[:, 1])
  
  if num_params_psd_per > 0:
    print "AP:", np.mean(sample[:, 3]/np.log(10.0)) + 2.0*np.log10(flux_scale), np.std(sample[:, 3]/np.log(10.0))
    print "nu0:", np.mean(sample[:, 4]/np.log(10.0)), np.std(sample[:, 4]/np.log(10.0))
    print "w:", np.mean(sample[:, 5]/np.log(10.0)), np.std(sample[:, 5]/np.log(10.0))

#=======================================================
# load observed light curve
#
#=======================================================
def load_lcdata():
  global lc, freqlc, psdlc
  global W, V, DT, time_media, flux_scale, flux_mean
  global parset
  
  lc = np.loadtxt(parset["FileName"])
  
  ts = lc[:, 0]
  fs  = copy.copy(lc[:, 1])
  if flag_endmatch:
    fs -= ((fs[-1] - fs[0])/(ts[-1] - ts[0]) * (ts - ts[0]) + fs[0])
  freqlc, psdlc = fft_psd(ts, fs)
  #freqlc, psdlc = np.loadtxt(os.path.dirname(parset["FileName"])+"/psd_"+os.path.basename(parset["FileName"]), unpack=True)
  
  cad = lc[1:, 0] - lc[0:-1, 0]
  cad_sort = np.sort(cad)
  time_cad_media = cad_sort[len(lc[:, 0])//2-1]
  DT = time_cad_media/W
  time_media = (lc[0, 0] + lc[-1, 0])/2.0
  flux_mean = np.mean(lc[:, 1])
  flux_scale = (np.max(lc[:, 1]) - np.min(lc[:, 1]))/2.0
  print "DT:", DT, time_media, flux_scale, flux_mean

#=======================================================
# load param  
#
#=======================================================
def load_param():
  global parset
  parset = {}
  fp = open("param", "r")
  for line in fp.readlines():
    line = line.lstrip()
    line = line.rstrip() + "  #"
    line = line.lstrip()
    if line[0] is not '#':
      arr = line.split()
      parset[arr[0]] = arr[1]
  
  parset['PSDType'] = int(parset['PSDType'])
  parset['V'] = float(parset['V'])
  parset['W'] = float(parset['W'])
  parset['FreqLimit'] = float(parset['FreqLimit'])
  parset['FlagEndMatch'] = int(parset['FlagEndMatch'])
  parset['FlagWhiteNoise'] = int(parset['FlagWhiteNoise'])
  parset['PeriodPSDProfType'] = int(parset['PeriodPSDProfType'])
  
  
  return

#=======================================================  
# calculate pb(TR)
#
#=======================================================
def pb_TR(doplot=False):
  global sample, freqlc, psdlc
  global lc, ncycle, tpmin
  global psd_period, psd_period_sqrt
  global ns
  
  TR = np.zeros(ns)
  TRobs = np.zeros(ns)
  
  tspan = lc[-1, 0] - lc[0, 0]
  idx = np.where((freqlc <= 1.0/tpmin) & (freqlc>= 1.0/(tspan/ncycle)))
  for k in range(ns):
    i = np.random.randint(sample.shape[0])
    ts, fs = genlc_psd_data(sample[i, :])

    if flag_endmatch:
      fs -= ((fs[-1] - fs[0])/(ts[-1] - ts[0]) * (ts - ts[0]) + fs[0])

    freq, psds = fft_psd(ts, fs)
    if num_params_psd_per == 0:
      psdtrue = psd_power_law(freq, sample[i, :num_params_psd_sto])*flux_scale**2
    else:
      psdtrue = (psd_power_law(freq, sample[i, :3]) + \
      psd_period(freq, sample[i, num_params_psd_sto:num_params_psd]))*flux_scale**2
      
    TR[k] = 2.0*np.max(psds[idx[0]]/psdtrue[idx[0]])
    TRobs[k] = 2.0*np.max(psdlc[idx[0]]/psdtrue[idx[0]])
    
    #plt.plot(freq, psds)
    #plt.plot(freq, psdtrue)
    #plt.plot(freqlc, psdlc, color='k')
    #plt.axvline(x=1.0/tpmin, ls='--')
    #plt.axvline(x=1.0/(tspan/ncycle),ls='--')    
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()
  
  pb = np.sum(TR>TRobs)*1.0/ns
  print "TR:", pb
  
  if doplot:
   plt.hist(TR, bins=50, normed=True, color='r')
   plt.hist(TRobs, bins=50, normed=True, color='b')
   plt.show()
   plt.close()

  return pb, TR, TRobs

#=======================================================
# calculate pb(TLS)
#
#=======================================================
def pb_TLS(doplot=False):
  global lc
  global ls_freq
  global ns
  
  ts = lc[:, 0]
  fs  = copy.copy(lc[:, 1])
  if flag_endmatch:
    fs -= ((fs[-1] - fs[0])/(ts[-1] - ts[0]) * (ts - ts[0]) + fs[0])
  lsp_data = LombScargle(ts, fs).power(ls_freq, normalization='standard')
  idxmax = np.argmax(lsp_data)
  
  period_obs = 1.0/ls_freq[idxmax]/365.0
  lspmax_obs = lsp_data[idxmax]
  print "LS data:", (period_obs, lspmax_obs)

  TLS_obs = lspmax_obs
  TLS = np.zeros(ns)
  
  for k in range(ns):
    i = np.random.randint(sample.shape[0])
    ts, fs = genlc_psd_data(sample[i, :])
    if flag_endmatch:
      fs -= ((fs[-1] - fs[0])/(ts[-1] - ts[0]) * (ts - ts[0]) + fs[0])
    lsp = LombScargle(ts, fs).power(ls_freq, normalization='standard')
    idxmax = np.argmax(lsp)
    period = 1.0/ls_freq[idxmax]/365.0
    TLS[k] = lsp[idxmax]

    #plt.plot(1.0/ls_freq, lsp)
    #plt.plot(1.0/ls_freq, lsp_data)
    #plt.show()
   
    
  pb_TLS = np.sum(TLS>TLS_obs)*1.0/ns
  print "TLS:", pb_TLS
  
  if doplot:
    plt.hist(TLS, bins=50)
    plt.axvline(x=TLS_obs, color='r')
    plt.show()
    plt.close()

  return pb_TLS, TLS, TLS_obs

#=======================================================
#calculate pb(TPDM)  
#
#=======================================================
def pb_TPDM(doplot=False):
  global lc, sample
  global pdm_scan
  global ns
  
  ts = lc[:, 0]
  fs  = copy.copy(lc[:, 1])

  if flag_endmatch:
    fs -= ((fs[-1] - fs[0])/(ts[-1] - ts[0]) * (ts - ts[0]) + fs[0])
  
  P = pyPDM.PyPDM(ts, fs)
  fpdm, tpdm_data = P.pdmEquiBinCover(5, 5, pdm_scan)
  idxmin = np.argmin(tpdm_data)
  print "PDM:", fpdm[idxmin], tpdm_data[idxmin]
  
  TPDM_obs = 1.0 - tpdm_data[idxmin]
  TPDM = np.zeros(ns)
  
  for k in range(ns):
    i = np.random.randint(sample.shape[0])
    ts, fs = genlc_psd_data(sample[i, :])
    fs -= ((fs[-1] - fs[0])/(ts[-1] - ts[0]) * (ts - ts[0]) + fs[0])
    P = pyPDM.PyPDM(ts, fs)
    fpdm, tpdm = P.pdmEquiBinCover(5, 5, pdm_scan)
    idxmin = np.argmin(tpdm)
    TPDM[k] = 1.0 - tpdm[idxmin]

    #plt.plot(fpdm, tpdm)
    #plt.plot(fpdm, tpdm_data)
    #plt.show()
  
  pb_TPDM = np.sum(TPDM>TPDM_obs)*1.0/ns
  print "TPDM:", pb_TPDM
  
  if doplot:
    plt.hist(TPDM, bins=50)
    plt.axvline(x=TPDM_obs, color='r')
    plt.show()
    plt.close()  
  return pb_TPDM, TPDM, TPDM_obs
    
  return 
  

def test_genlc():
  ts, fs = genlc(sample[0, :])
  plt.plot(ts, fs)

  recon = np.loadtxt("data/recon_all.txt")
  offset = recon.shape[0]/sample.shape[0]
  plt.plot(recon[:offset, 0], recon[:offset, 1])

  plt.show()
  
if __name__=="__main__":
  global sample, lc, freqlc, psdlc
  global freq_limit, num_params_psd, num_params_psd_sto, num_params_psd_per, nd_sim
  global W, V, DT, time_media, flux_scale, flux_mean
  global parset
  global pdm_scan, ls_freq, ncycle, tpmin
  global psd_period, psd_period_sqrt
  global flag_endmatch
  global ns
  
  ns = 1000
  ncycle = 3.
  tpmin = 100.0
  flag_endmatch = True
  print flag_endmatch

  load_param()
  
  if parset['PeriodPSDProfType'] == 0:
    psd_period = psd_period_gaussian
    psd_period_sqrt = psd_period_sqrt_gaussian
  else:
    psd_period = psd_period_lorentz
    psd_period_sqrt = psd_period_sqrt_lorentz

  num_params_psd_sto = 3
  if(parset['PSDType'] >= 3):
    num_params_psd_per = 3
  else:
    num_params_psd_per = 0
    
  num_params_psd = (num_params_psd_sto+num_params_psd_per) 
  freq_limit = parset["FreqLimit"]
  V = parset['V']
  W = parset['W']
  
  load_lcdata()
  load_sample()
  
  tspan = lc[-1, 0] - lc[0, 0]
  print "Tspan:", tspan/365.0
  ls_freq = np.logspace(np.log10(1.0/tpmin), np.log10(1.0/(tspan/ncycle)), 500)
  pdm_scan = pyPDM.Scanner(minVal=tpmin, maxVal=tspan/ncycle, dVal=(tspan/ncycle-tpmin)/500, mode="period")
  
  if num_params_psd_per == 0:
    nd_sim = sample.shape[1] - num_params_psd
  else:
    nd_sim = ((sample.shape[1] - num_params_psd) * 2)//3
  
  print sample.shape[1], nd_sim, num_params_psd

  #test_genlc()

  pbTR, TR, TR_obs = pb_TR(doplot=True)
  pbTLS, TLS, TLS_obs = pb_TLS(doplot=True)
  pbTPDM, TPDM, TPDM_obs = pb_TPDM(doplot=True)

  np.savetxt("TR.txt", np.stack((TR, TR_obs), axis=-1), fmt="%15.5f  %15.5f")
  np.savetxt("TLS.txt", np.hstack((TLS, TLS_obs)), fmt="%15.5f")
  np.savetxt("TPDM.txt", np.hstack((TPDM, TPDM_obs)), fmt="%15.5f")


