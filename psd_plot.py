#
# a python script for showing results.
#
# Yan-Rong Li, liyanrong@mail.ihep.ac.cn
# Nov 3, 2018

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import configparser as cfgpars 

def read_params(fname):
  """
  read parameter file
  """
  config = cfgpars.RawConfigParser(delimiters=' ', comment_prefixes='#', inline_comment_prefixes='#', \
      default_section=cfgpars.DEFAULTSECT, empty_lines_in_values=False)
  
  with open(fname) as f:
      file_content = '[dump]\n' + f.read()
  
  config.read_string(file_content)
  params = config["dump"]
  
  return params

def psd_harmonic(params, freq, arg):
  """
  PSD function for harmonic.
  """
  harmonic_term_num = int(params["PSDType"])
  psd = np.zeros(freq.shape)
  S1 = np.exp(arg[0])
  omega1_sqr = np.exp(2.0*arg[1])
  psd[:] += S1 * omega1_sqr**2/((freq**2 - omega1_sqr)**2 + 2*omega1_sqr * freq**2)
  noise = np.exp(arg[2 + (harmonic_term_num-1) * 3])

  for j in range(harmonic_term_num-1):
    S2 = np.exp(arg[2+j*3])
    omega2_sqr = np.exp(2.0*arg[2+j*3+1])
    Q_sqr = (np.exp(arg[2+j*3+2])-1.0)**2/4.0
    psd[:] += S2 * omega2_sqr**2/((freq**2 - omega2_sqr)**2 + omega2_sqr * freq**2/Q_sqr)
    
  psd[:] *= np.sqrt(2.0/np.pi)
   
  return psd + noise

def harmonic_plot(params):
  """
  plot results for harmonic model
  """
  lc = np.loadtxt(params["FileName"])
  lc_recon = np.loadtxt("data/recon_mean.txt")
  harmonic_term_num = int(params["PSDType"])
  num_params = 2 + 3*(harmonic_term_num - 1) + 1
  sample = np.loadtxt("data/posterior_sample.txt", usecols=(range(num_params)))
  scale = (np.max(lc[:, 1]) - np.min(lc[:, 1]))/2.0
  psd_data = np.loadtxt(os.path.dirname(params["FileName"])+"/psd_"+os.path.basename(params["FileName"]))

  freq = np.logspace(np.log10(1.0/(lc[-1, 0]-lc[0, 0])), np.log10(len(lc[:, 0])/(lc[-1, 0]-lc[0, 0])), 1000)
  psd = np.zeros((sample.shape[0], 1000))
  
  fig = plt.figure(figsize=(15, 4))
  ax = fig.add_subplot(121)
  
  ax.errorbar(lc[:, 0], lc[:, 1], yerr = lc[:, 2], ls='none', marker='o', color='k')
  ax.plot(lc_recon[:, 0], lc_recon[:, 1])
  ax.fill_between(lc_recon[:, 0], y1 = lc_recon[:, 1] - lc_recon[:, 2], y2 = lc_recon[:, 1] + lc_recon[:, 2], zorder=0)
  ax.set_xlim(lc[0, 0]-0.1*(lc[-1, 0]-lc[0, 0]), lc[-1, 0]+0.1*(lc[-1, 0]-lc[0, 0]))
  fmax = np.max(lc[:, 1])
  fmin = np.min(lc[:, 1])
  ax.set_ylim(fmin - 0.1*(fmax-fmin), fmax + 0.1*(fmax-fmin))
  ax.set_xlabel('Time')
  ax.set_ylabel('Flux')
  
  ax = fig.add_subplot(122)
  for i in range(sample.shape[0]):
    psd[i, :] = psd_harmonic(params, freq, sample[i, :])
  
  psd_mean = np.mean(np.log10(psd), axis=0)
  psd_upper, psd_lower = np.percentile(np.log10(psd), [(100.0-68.3)/2.0, 100.0-(100.0-68.3)/2.0], axis=0)
  psd_best = psd_harmonic(params, freq, np.mean(sample, axis=0))
  
  ax.fill_between(freq, y1 = 10.0**psd_upper*scale**2, y2 = 10.0**psd_lower*scale**2, alpha=0.5)
  ax.plot(freq, 10.0**psd_mean*scale**2, color='b')
  ax.plot(psd_data[:, 0], psd_data[:, 1])
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlabel('Frequency')
  ax.set_ylabel('PSD')
  plt.show()

def psd_carma(params, freq, arg):
  """
  psd function for carma process
  """
  pstr = params["PSDType"].split(":")
  carma_p = int(pstr[0])
  carma_q = int(pstr[1])
  num_params = carma_p + carma_q + 1 + 1
  
  roots = np.zeros(carma_p, dtype=complex)
  for i in range(carma_p//2):
    a = np.exp(arg[1+i*2])
    b = np.exp(arg[1+i*2+1])
    roots[i*2] = -2.0*np.pi*(a + b*1j)
    roots[i*2+1] = -2.0*np.pi*(a - b*1j)
  if carma_p%2 == 1:
    roots[carma_p-1] = -2.0*np.pi* np.exp(arg[1+carma_p -1])
  
  coefs = np.poly(roots)
  ar_coefs=np.zeros(carma_p+1)
  ar_coefs[:] = coefs[::-1].real
  
  ma_coefs = np.zeros(carma_p+1)
  ma_coefs[0] = 1.0
  for i in range(1, carma_q+1):
    ma_coefs[i] = np.exp(arg[1+carma_p + i -1])
  for i in range(carma_q+1, carma_p+1):
    ma_coefs[i] = 0.0
  
  ar_poly = np.zeros(len(freq), dtype=complex)
  ma_poly = np.zeros(len(freq), dtype=complex)
  for k in range(carma_p+1):
    tmp = (2.0*np.pi*freq*1j)**(k)
    ar_poly += ar_coefs[k] * tmp
    ma_poly += ma_coefs[k] * tmp
  
  psd = np.exp(2.0*arg[0]) * np.abs(ma_poly)**2/np.abs(ar_poly)**2 

  return psd 
  

def carma_plot(params):
  """
  plot results for carma model
  """
  pstr = params["PSDType"].split(":")
  carma_p = int(pstr[0])
  carma_q = int(pstr[1])
  num_params = carma_p + carma_q + 1 + 1
  print(num_params)

  lc = np.loadtxt(params["FileName"])
  lc_recon = np.loadtxt("data/recon_mean.txt")
  sample = np.loadtxt("data/posterior_sample.txt", usecols=(range(num_params)))
  
  scale = (np.max(lc[:, 1]) - np.min(lc[:, 1]))/2.0
  psd_data = np.loadtxt(os.path.dirname(params["FileName"])+"/psd_"+os.path.basename(params["FileName"]))
  
  freq = np.logspace(np.log10(1.0/(lc[-1, 0]-lc[0, 0])), np.log10(len(lc[:, 0])/(lc[-1, 0]-lc[0, 0])), 1000)
  psd = np.zeros((sample.shape[0], 1000))
  
  fig = plt.figure(figsize=(15, 4))
  ax = fig.add_subplot(121)
  
  ax.errorbar(lc[:, 0], lc[:, 1], yerr = lc[:, 2], ls='none', marker='o', color='k')
  ax.plot(lc_recon[:, 0], lc_recon[:, 1])
  ax.fill_between(lc_recon[:, 0], y1 = lc_recon[:, 1] - lc_recon[:, 2], y2 = lc_recon[:, 1] + lc_recon[:, 2], zorder=0)
  ax.set_xlim(lc[0, 0]-0.1*(lc[-1, 0]-lc[0, 0]), lc[-1, 0]+0.1*(lc[-1, 0]-lc[0, 0]))
  fmax = np.max(lc[:, 1])
  fmin = np.min(lc[:, 1])
  ax.set_ylim(fmin - 0.1*(fmax-fmin), fmax + 0.1*(fmax-fmin))
  ax.set_xlabel('Time')
  ax.set_ylabel('Flux')
  
  ax = fig.add_subplot(122)
  for i in range(sample.shape[0]):
    psd[i, :] = psd_carma(params, freq, sample[i, :])
  
  psd_mean = np.mean(np.log10(psd), axis=0)
  psd_upper, psd_lower = np.percentile(np.log10(psd), [(100.0-68.3)/2.0, 100.0-(100.0-68.3)/2.0], axis=0)
  psd_best = psd_carma(params, freq, np.mean(sample, axis=0))
  
  ax.fill_between(freq, y1 = 10.0**psd_upper*scale**2, y2 = 10.0**psd_lower*scale**2, alpha=0.5)
  ax.plot(freq, 10.0**psd_mean*scale**2, color='b')
  ax.plot(psd_data[:, 0], psd_data[:, 1]*scale**2)
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlabel('Frequency')
  ax.set_ylabel('PSD')
  plt.show()

def psd_simple(params, freq, arg):

  num_params_psd = 0
  if params["PSDType"] == "0":
    A = np.exp(arg[0])
    alpha = arg[1]
    noise = np.exp(arg[2])
    psd = A * (freq)**(-alpha) + noise
    num_params_psd = 3
  elif params["PSDType"] == "1":
    A, fknee, noise = np.exp(arg[0:3])
    psd = A / (1.0 + (freq/fknee)**2 ) + noise
    num_params_psd = 3
  else:
    A = np.exp(arg[0])
    alpha_hi = arg[1]
    alpha_lo = arg[1] - arg[2]
    fc, noise = np.exp(arg[3:5])
    psd = np.zeros(len(freq))
    idx = (freq >= fc)
    psd[idx] = A * (freq[idx]/fc)**(-alpha_hi) + noise
    idx = (freq < fc)
    psd[idx] = A * (freq[idx]/fc)**(-alpha_lo) + noise
    num_params_psd = 5
  
  if params["PSDPeriodModel"].lower() == "delta":
    Ap = np.exp(arg[num_params_psd])
    nu0 = np.exp(arg[num_params_psd+1])
    phi = np.exp(arg[num_params_psd+2])
    i0 = int( (np.log10(nu0) - np.log10(freq[0]))/(np.log10(freq[1]) - np.log10(freq[0])))
    psd[i0] += Ap
  elif params["PSDPeriodModel"].lower() == "gaussian":
    Ap, center, sig = np.exp(arg[num_params_psd:num_params_psd+3])
    psd += Ap/np.sqrt(2.0*np.pi)/sig * np.exp(-0.5*(freq - center)**2/sig**2)
  elif params["PSDPeriodModel"].lower() == "lorentzian":
    Ap, center, sig = np.exp(arg[num_params_psd:num_params_psd+3])
    psd += Ap/np.pi * sig/(sig*sig + (freq - center)**2)

  return psd

def simple_plot(params):
  lc = np.loadtxt(params["FileName"])
  lc_recon = np.loadtxt("data/recon_mean.txt")
  psdtype = int(params["PSDType"])
  if params["PSDType"] == "0":
    num_params = 3
  elif params["PSDType"] == "1":
    num_params = 3
  else:
    num_params = 5

  if params["PSDPeriodModel"].lower() == "delta":
    num_params += 3
  elif params["PSDPeriodModel"].lower() == "gaussian":
    num_params += 3
  elif params["PSDPeriodModel"].lower() == "lorentzian":
    num_params += 3
  else:
    num_params += 0

  sample = np.loadtxt("data/posterior_sample.txt", usecols=(range(num_params)))

  scale = (np.max(lc[:, 1]) - np.min(lc[:, 1]))/2.0
  psd_data = np.loadtxt(os.path.dirname(params["FileName"])+"/psd_"+os.path.basename(params["FileName"]))

  freq = np.logspace(np.log10(1.0/(lc[-1, 0]-lc[0, 0])), np.log10(len(lc[:, 0])/(lc[-1, 0]-lc[0, 0])), 1000)
  psd = np.zeros((sample.shape[0], 1000))
  
  fig = plt.figure(figsize=(15, 4))
  ax = fig.add_subplot(121)
  
  ax.errorbar(lc[:, 0], lc[:, 1], yerr = lc[:, 2], ls='none', marker='o', color='k')
  ax.plot(lc_recon[:, 0], lc_recon[:, 1])
  ax.fill_between(lc_recon[:, 0], y1 = lc_recon[:, 1] - lc_recon[:, 2], y2 = lc_recon[:, 1] + lc_recon[:, 2], zorder=0)
  ax.set_xlim(lc[0, 0]-0.1*(lc[-1, 0]-lc[0, 0]), lc[-1, 0]+0.1*(lc[-1, 0]-lc[0, 0]))
  fmax = np.max(lc[:, 1])
  fmin = np.min(lc[:, 1])
  ax.set_ylim(fmin - 0.1*(fmax-fmin), fmax + 0.1*(fmax-fmin))
  ax.set_xlabel('Time')
  ax.set_ylabel('Flux')
  
  ax = fig.add_subplot(122)
  for i in range(sample.shape[0]):
    psd[i, :] = psd_simple(params, freq, sample[i, :])
  
  psd_mean = np.mean(np.log10(psd), axis=0)
  psd_upper, psd_lower = np.percentile(np.log10(psd), [(100.0-68.3)/2.0, 100.0-(100.0-68.3)/2.0], axis=0)
  psd_best = psd_harmonic(params, freq, np.mean(sample, axis=0))
  
  ax.fill_between(freq, y1 = 10.0**psd_upper*scale**2, y2 = 10.0**psd_lower*scale**2, alpha=0.5)
  ax.plot(freq, 10.0**psd_mean*scale**2, color='b')
  ax.plot(psd_data[:, 0], psd_data[:, 1])
  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.set_xlabel('Frequency')
  ax.set_ylabel('PSD')
  plt.show()

def do_plot(params):
  if params["PSDModel"].lower() == "harmonic":
    harmonic_plot(params)
  elif params["PSDModel"].lower() == "carma":
    carma_plot(params)
  elif params["PSDModel"].lower() == "simple":
    simple_plot(params)
    

if __name__=="__main__":
  
  if len(sys.argv) < 2:
    raise Exception("Please specify a param file as: 'python psd_plot.py param'.")

  params = read_params(sys.argv[1])
  do_plot(params)
  


