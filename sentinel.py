#!/usr/bin/env python

'''
RSE Lewis et al. 2012 experiment #1
##############################################
This code generates a synthetic set of Sentinel-2 MSI
data and tests the DA system using them
In this first experiment, we generate a clean (no noise)
dataset with full bandpass sampling and attempt to recover
the model parameters. This uses no regularisation, but has
full temporal coverage.
'''
import numpy as np
import ConfigParser
import os
import errno
try:
  import ephem   # easy_install pyephem
except:
  print 'Error loading ephem package'
  print 'Try installing it'
  print 'e.g. easy_install pyephem'

import datetime
import eoldas
import tempfile
import pylab as plt


class Sentinel():
  '''
   This code generates a synthetic set of Sentinel-2 MSI
   data following Lewis et al. 2012 (RSE) but has some options.

   e.g. use:

   import sentinel 
   import numpy as np

   s = sentinel.Sentinel()
   # generate parameter dataset
   s.parameters(np.arange(1,366),'input/truth.dat')
   # generate noise free observations
   s.fwdModel('input/truth.dat','input/sentinelClean.dat')
   # add noise to the observations
   s.addNoiseToObservations('input/sentinelClean.dat','input/sentinel.dat')
   # solve for parameter estimate, using one date at a time
   # e.g. test for clean data starting from the correct value
   # this should produce the correct values
   s.solveSingle('input/sentinelClean.dat','output/sentinelClean1.dat',initial='input/truth.dat')
   truth = s.loadData('input/truth.dat')
   est   = s.loadData('output/sentinelClean1.dat')
 
   # e.g. test for clean data not starting from the correct value
   s.solveSingle('input/sentinelClean.dat','output/sentinelClean2.dat')
   # e.g. test normal
   s.solveSingle('input/sentinel.dat','output/sentinel.dat')

  '''
  def __init__(self):
    # improve this later 
    # set some default parameter values
    self.gamma = 1
    self.xlai = np.exp(-3.0/2.)
    self.xhc = 5
    self.rpl = 0.01
    self.xkab = np.exp(-50/100.)
    self.scenesc = 0.0
    self.xkw = np.exp(-0.01*50.)
    self.xkm = np.exp(-100.*0.01)
    self.xleafn = 1.5
    self.xs1 = 1.0
    self.xs2 = 0
    self.xs3 = 0
    self.xs4 = 0
    self.lad = 5
    self.vary = {'xlai':True,'xkab':True,'xkw':True,'xkm':True,'xs1':True}

  def loadData(self,file):
    '''
    Load data from file into a convenient format
    using keys from the header in the file
    '''
    theader,tdata = self.readParameters(file)
    out = {}
    for i,k in enumerate(theader.split()):
      out[k] = tdata[:,i] 
    return out

  def parameters(self,doys,truthFile,confFile=None):
    '''
    Generate a 'parameter' dataset into the file self.truthFile
    ('input/truth.dat' by default) based on the temporal functions 
    in Lewis et al. (2012).

    Default parameter values are picked up from
      self.gamma,self.xlai,self.xhc,self.rpl,self.xkab,self.scenesc,self.xkw
      self.xkm,self.xleafn,self.xs1,self.xs2,self.xs3,self.xs4,self.lad
    which have default values set upon initialisation.

    doys        : a np array of day numbers

    Options:
      confFile  : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.

    '''
    try:
      self.doys = doys
    except:
       self.error('failure setting doys field in self.parameters()')
       return

    self.truthFile = truthFile

    # read the config file to get the defaults
    config = ConfigParser.RawConfigParser()
    if confFile != None:
      self.confFile = confFile
    try:
      config.read(self.confFile)
    except:
      self.writeConf(confFile=confFile)
      config.read(self.confFile)

    self.names = [i.strip() for i in config.get('parameter','names').split(',')]
    self.params = {}

    cmd = '--parameter.x.default=%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d' \
          %(self.gamma,self.xlai,self.xhc,self.rpl,self.xkab,self.scenesc,\
            self.xkw,self.xkm,self.xleafn,self.xs1,self.xs2,self.xs3,self.xs4,self.lad)

    self.paramlist = cmd.split('=')[1].split(',')

    for (j,i) in enumerate(self.names):
      self.params[i] = float(self.paramlist[j])

    self.data = np.zeros([len(doys),len(self.params)+2])

    # set data array with all the values to write out
    self.data[:,0] = self.doys
    self.data[:,1] = 1
    self.datastr = 'time mask'
    for (n,i) in enumerate(self.names):
      t = (self.doys-1)/365.
      self.datastr = self.datastr + ' %s'%i
      if i == 'gamma_time':
        self.data[:,n+2] = self.gamma
      elif self.vary['xlai']  and i == 'xlai':
        self.data[:,n+2] = 0.21 + 3.51 * (np.sin(np.pi*t)**5)
        self.data[:,n+2] = np.exp(-self.data[:,n+2]/2.)
      elif self.vary['xkab']  and i == 'xkab':
        w = np.where(t<=0.5)[0]
        self.data[w,n+2] = 10.5 + 208.7*t[w]
        w = np.where(t>0.5)[0]
        self.data[w,n+2] = 219.2 - 208.7*t[w]
        self.data[:,n+2] = np.exp(-self.data[:,n+2]/100.)
      elif self.vary['xkw']  and i == 'xkw':
	#data[:,n+2] = 0.0068 + 0.002*np.sin(np.pi * t+0.1) * np.sin(6*np.pi*t + 0.1)
        # inconsistent in the paper use:
        self.data[:,n+2] = 0.068/5 + 0.01*np.sin(np.pi * t+0.1) * np.sin(6*np.pi*t + 0.1)
        self.data[:,n+2] = np.exp(-self.data[:,n+2]*50.)
      elif self.vary['xs1']  and i == 'xs1':
        # difft soil model here ...  so scale
	self.data[:,n+2] = 2.5*(0.2 + 0.18*np.sin(np.pi*t) * np.sin(6*np.pi*t))
        #data[:,n+2] = 2.5*(0.18*np.sin(np.pi*t) * np.sin(6*np.pi*t))
      else:
        self.data[:,n+2] = self.params[i]

    self.writeParameters(self.truthFile,self.datastr,self.data)
 
  def getBands(self,confFile=None):
    '''
    Read a configuration file and set self.bands, self.control, self.mask

    Options:
      confFile   : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.

    '''
    # read the config file to get the defaults
    config = ConfigParser.RawConfigParser()
    if confFile != None:
      self.confFile = confFile
    try:
      config.read(self.confFile)
    except:
      self.writeConf(confFile=confFile)
      config.read(self.confFile)

    # wavebands
    self.bands = [i.strip().replace("'","").replace("[","").replace("]","") \
                  for i in config.get('operator.obs.y','names').split(',')]
    self.control = eval(config.get('operator.obs.y','control'))
 
  def fwdModel(self,ifile,ofile, \
                    lat='50:0', lon='0:0', year='2011', maxVza=15.0, \
                    minSD=0.004, maxSD=0.01,fullBand=False, \
                    confFile=None,every=5,prop=1.0,WINDOW=1):
    '''
    Generate a fwd modelling (reflectrance data) using the parameters in ifile
    writing out to ofile.

    ifile    : input data (parameters) file
    ofile    : output reflectance data file

    Options:
      lat        : latitude (default '50:0') (see ephem)
      lon        : longitude (default ''0:0') (see ephem)
      year       : as int or string (default '2011')
      maxVza     : maximum view zenith angle assumed (15 degrees default)
      minSD      : minimum noise
      maxSD      : maximum noise. The uncertainty in each waveband is scaled linearly
                   with wavelength between minSD and maxSD 
      fullBand   : set True if you want full band pass else False (default False)
                   Note that its much slower to set this True.
      confFile   : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.
      every      : sample output every 'every' doys (default 5)
      prop       : proportion of clear sample days (default 1.0)
      WINDOW     : size of smoothing kernel to induce temporal correlation in 
                   data gaps if prop < 1
    '''
    this = self.loadData(ifile)
    self.doys = this['time']

    # read the config file to get the defaults
    self.getBands(confFile=confFile)

    self.mask = np.zeros_like(self.doys).astype(bool)
    self.datastr = 'time '

    # sub sample
    self.sdoys = self.doys[self.doys%every == 1]

    # apply gaps
    WINDOW = 1
    weightings = np.repeat(1.0, WINDOW) / WINDOW
    l = len(self.sdoys)
    xx = np.convolve(np.random.rand(l*100),weightings,'valid')[WINDOW:WINDOW+l]

    maxx = sorted(xx)[:int(len(xx)*prop)]
    mask = np.in1d(xx,maxx)
    self.sdoys = self.sdoys[mask]

    # generate data
    l = len(self.sdoys)
    self.sdata = np.zeros([l,1 + len(self.bands)*2+len(self.control)])
    self.sdata[:,0] = self.sdoys
    o = ephem.Observer()
    o.lat, o.long, o.date = lat, lon , datetime.datetime(int(year), 1, 1, 10, 30)
    sun = ephem.Sun(o)
    dd = o.date

    for (n,i) in enumerate(self.control):
      self.datastr += ' %s'%i
      if i == 'mask':
        self.sdata[:,n+1] = 1
      elif i == 'vaa' or i == 'saa':
        self.sdata[:,n+1] = np.random.rand(l)*360.
      elif i == 'vza':
        self.sdata[:,n+1] = np.random.rand(l)*maxVza
      else:
        dates = dd + self.sdoys
        for (m,j) in enumerate(dates):
          o.date = j
          sun = ephem.Sun(o)
          self.sdata[m,n+1] = 90 - float(sun.alt) * 180./np.pi

    for (n,i) in enumerate(self.bands):
      self.datastr += ' %s'%i

    i = self.bands[0]
    this0 = this = np.array([float(j) for j in i.split('-')]).mean()
    i = self.bands[-1]
    this1 = np.array([float(j) for j in i.split('-')]).mean()

    ostr = ''
    nstr = ''
    # sort the sd info
    for (n,i) in enumerate(self.bands):
      ostr = ostr + "%12.6f"%0.0
      this = minSD + (maxSD - minSD) * (np.array([float(j) for j in i.split('-')]).mean() - this0)/(this1 - this0)
      nstr = nstr + "%12.6f"%this
      self.datastr += ' sd-%s'%i
      self.sdata[:,1+len(self.control)+len(self.bands)+ n] = this

    # write out a first version of the observations file   
    self.writeParameters(ofile,self.datastr,self.sdata) 

    ############################################
    # Generate the clean synthetic observations
    ############################################

    # we ingest the state vector data (ifile)
    # and the angle / observation data (ofile)
    # and produce fwdFile that has the (clean) fwd observations
    # but doesn't have the uncertainty information
    # Note that we can switch use_median False, which
    # means that the full bandpass is used in fwd modelling
    fwdFile = ofile + '.tmp'
    fwdData = ifile + '.tmp'
    if fullBand:
      option = "--operator.obs.rt_model.no_use_median "
    else:
      option = ""

    fwdCmd = "eoldas --conf=config_files/eoldas_config.conf --conf=config_files/sentinel0.conf  --passer " + \
                "--operator.obs.y.result.filename=%s "%fwdFile + \
                "--operator.obs.y.state=%s "%ofile + \
                "--logfile=logs/rseFwd.log %s "%option + \
                "--parameter.result.filename=%s "%fwdData + \
                "--conf=config_files/sentinel1.conf --parameter.x.state=%s "%ifile

    self = eoldas.eoldas(fwdCmd)
    # write fwd file
    self.solver.write()
    self.solver.writeHx()

    # we now put the observation uncertainty information in
    # combination with the observation data in fwdFile
    # to give the clean observations in op
    ip = fwdFile
    op = ofile
    f = open(ip,'r')
    f2 = open(op,'w')

    this = f.readlines()
    that = [i.replace(ostr,nstr)  for i in this]

    f2.writelines(that)
    f.close()
    f2.close()


  def solveRegular(self,ifile,ofile,modelOrder=1,gamma=None,initial=None,confFile=None):
    '''
    Use eoldas to solve for parameter estimates using 
    observations in ifile, writing result to ofile.

    In this case, there is no temporal constraint, so only a very weak
    prior (see configuration file) is used and we solve for each day's observation 
    sequentially.

    ifile        : input observations file
    ofile        : output parameter file

    Options:
      modelOrder : regularisation model order (e.g. 1 or 2 (default 1)
      gamma      : gamma value (default whatever set on initialisation or 
                   in initial file)
      initial    : starting point (default is None, so start at the priors
                   defined in the configuration file)

      confFile   : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.

    '''
    self.getBands(confFile=confFile)
    if gamma:
      self.gamma = gamma

    cmd = '--parameter.x.default=%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d' \
          %(self.gamma,self.xlai,self.xhc,self.rpl,self.xkab,self.scenesc,\
            self.xkw,self.xkm,self.xleafn,self.xs1,self.xs2,self.xs3,self.xs4,self.lad)

    self.paramlist = cmd.split('=')[1].split(',')

    #set up a new conf to give this
    str = '''
[operator]
modelt.name=DModel_Operator
modelt.datatypes = x

[operator.modelt.x]
names = $parameter.names
sd = [1.0]*len($operator.modelt.x.names)
datatype = x

[operator.modelt.rt_model]
model_order=%s
wraparound=periodic,365

    '''%(modelOrder)
    
    conf2 = 'tmp/c2'
    self.mkdir(conf2)
    cmd += ' --conf=%s '%conf2
    open(conf2,'w').write(str)

    if initial:
      conf2 = 'tmp/c1'
      initfile = 'tmp/c3'
      self.mkdir(initfile)
      #set up a new conf to give this
      str = '''
[parameter.x]
state = %s
    '''%(initfile)      
      open(conf2,'w').write(str)
      iheader,idata = self.readParameters(initial)
      isdata = {}
      for i,h in enumerate(iheader.split()):
        isdata[h] = idata[:,i]
      idoys = isdata['time']
      cmd += ' --conf=%s '%conf2
      # fix the gamma value in this file to what we want it to be
      if gamma:
        isdata['gamma_time'] = self.gamma
      for i,h in enumerate(iheader.split()):
        idata[:,i] = isdata[h]
      self.writeParameters(initfile,iheader,idata)

    # read the input file
    header,data = self.readParameters(ifile)
    sdata = {}
    for i,h in enumerate(header.split()):
      sdata[h] = data[:,i]
    doys = sdata['time']

    opfile = ofile + '_result.dat'
    fwdFile = ofile + '_fwd.dat'
    pFile = opfile + '_prior.dat'

    fwdCmd = "eoldas --conf=config_files/eoldas_config.conf --conf=%s %s "%(self.confFile,cmd) + \
          "--operator.obs.y.result.filename=%s "%fwdFile + \
          "--operator.obs.y.state=%s "%ifile + \
          "--operator.prior.y.result=%s "%pFile + \
          "--logfile=logs/rseSolve.log --no_init_test --plotmod=50 --optimisation.gtol=1e-10 " + \
          "--parameter.result.filename=%s "%opfile

    eo = eoldas.eoldas(fwdCmd)
    eo.solve(write=True,unc=True)
                                      

  def solveSingle(self,ifile,ofile,initial=None,confFile=None):
    '''
    Use eoldas to solve for parameter estimates using 
    observations in ifile, writing result to ofile.

    In this case, there is no temporal constraint, so only a very weak
    prior (see configuration file) is used and we solve for each day's observation 
    sequentially.

    ifile        : input observations file
    ofile        : output parameter file

    Options:
      initial    : starting point (default is None, so start at the priors
                   defined in the configuration file)

      confFile   : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.

    '''
    self.getBands(confFile=confFile)

    cmd = '--parameter.x.default=%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d' \
          %(self.gamma,self.xlai,self.xhc,self.rpl,self.xkab,self.scenesc,\
            self.xkw,self.xkm,self.xleafn,self.xs1,self.xs2,self.xs3,self.xs4,self.lad)

    self.paramlist = cmd.split('=')[1].split(',')

    if initial:
      #set up a new conf to give this
      str = '''
[parameter.x]
state = %s
      '''%initial
      conf2 = 'tmp/c4'
      self.mkdir(conf2)
      open(conf2,'w').write(str)
      iheader,idata = self.readParameters(initial)
      isdata = {}
      for i,h in enumerate(iheader.split()):
        isdata[h] = idata[:,i]
      idoys = isdata['time']
      cmd += ' --conf=%s '%conf2
 
    # read the input file
    header,data = self.readParameters(ifile)
    sdata = {}
    for i,h in enumerate(header.split()):
      sdata[h] = data[:,i]
    doys = sdata['time']

    costs = []
    # something wrong here I think, but doesnt matter as we dont want them anyway
    # so effectively sends them to /dev/null I think
    opfile = tempfile.NamedTemporaryFile().name
    fwdFile = tempfile.NamedTemporaryFile().name
    pFile = tempfile.NamedTemporaryFile().name

    for i,d in enumerate(doys):
      if i == 0:
        fwdCmd = "eoldas --conf=config_files/eoldas_config.conf --conf=%s %s "%(self.confFile,cmd) + \
          "--operator.obs.y.result.filename=%s "%fwdFile + \
          "--operator.obs.y.state=%s "%ifile + \
          "--operator.prior.y.result=%s "%pFile + \
          "--logfile=logs/rseSolveSingle.log --no_init_test --plotmod=1e20 --no_doplot --optimisation.gtol=1e-10 " + \
          "--parameter.limits=[[%d,%d,1]] "%(d,d) + \
          "--parameter.result.filename=%s "%opfile

        eo = eoldas.eoldas(fwdCmd)
        prior = truth = eo.solver.root.x.state
        #eo.solver.root.x.state = prior
        eo.solve(write=True,unc=True)
        # get  the tmp op file
        theader,tdata = self.readParameters(opfile)
        data = np.zeros((len(doys),len(theader.split()))) 
        data[i,:] = tdata
        this = np.array([eo.solver.root.x.state[0],eo.solver.root.x.sd]).flatten()
        npm = len(this)
      else:
        # just load the new observation
        # into eo.solver.root.operators[1].y.state
        for k in xrange(len(eo.solver.root.operators)):
          if eo.solver.root.operators[k].thisname == 'eoldas.solver.eoldas.solver-obs':
            for j,b in enumerate(self.bands):
              eo.solver.root.operators[k].y.state[0][j] = sdata[b].flatten()[i]
            # load angles & mask (control data)
            # need to indicate a reload needed for metadata
            eo.solver.root.operators[k].isLoaded = False
            for j,b in enumerate(eo.solver.root.operators[k].y_meta.control):
              eo.solver.root.operators[k].y.control[0][j] = sdata[b][i]
        if initial == None:
          eo.solver.root.x.state[0,1:] = 0.5*(prior[0,1:] + eo.solver.root.x.state[0,1:])
        else:
          thisd = np.where(np.in1d(idoys,d))[0]
          for j,nn in enumerate(eo.solver.names):
            eo.solver.root.x.state[0,j] = isdata[nn][thisd][0]
        truth = eo.solver.root.x.state[0].copy()
        eo.solve(write=False,unc=True)
        data[i,-npm:] = np.array([eo.solver.root.x.state[0],eo.solver.root.x.sd]).flatten()
      try:
        costs.append(eo.solver.min_cost[0])
      except:
        costs.append(eo.solver.min_cost)
      #print i,d,'cost',eo.solver.min_cost[0]
      #print 'data ',eo.solver.root.x.state[0]
      if initial:
        print 'truth',truth
    data[:,0] = doys
    self.data = data
    self.header = header
    self.eoldas = eo
    self.costs = costs
    self.writeParameters(ofile,theader,data)

  def plotParameters(self,truth,inferred,confFile=None):
    '''
    Plot truth and inferred parameters

      confFile   : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.

    '''
    self.getBands(confFile=confFile)


  def addNoiseToObservations(self,ifile,ofile,nMag=1.0,confFile=None):
    '''
    Add noise to observations in PARAMETERS file ifile
    and write to ofile

    Options:
      nMag    : magnify the noise by nMag (default 1.0)
      confFile   : configuration file (default config_files/sentinel0.conf).
                  If this file doesnt exist, it will be generated from 
                  self.confTxt. If that doesn't exist, self.generateConfTxt()
                  is invoked to provide a default.
    '''
    # read the config file to get the defaults
    self.getBands(confFile=confFile)

    header,data = self.readParameters(ifile)
    sdata = {}
    for i,h in enumerate(header.split()):
      sdata[h] = data[:,i]

    for b in self.bands:
      sdata['sd-' + b] *= nMag
      sdata[b] = np.random.normal(sdata[b],sdata['sd-' + b])

    for i,h in enumerate(header.split()):
      data[:,i] = sdata[h]
    self.writeParameters(ofile,header,data) 

  def readParameters(self,file):
    '''
    Read a parameters file

    return header and data
    '''
    header = open(file).readline().replace('#PARAMETERS ','')
    data = np.loadtxt(file)
    return header,data

  def writeParameters(self,file,header,data):
    '''
    Write a PARAMETERS file
    
    file    : filename
    header  : the data field names (as one string)
    data    : the data
    '''
    # make sure directory exists & get rid of any newline
    header = header.replace('\n','')
    self.mkdir(file)
    f = open(file,'w')
    f.write('#PARAMETERS %s\n'%header)

    ff2 = tempfile.NamedTemporaryFile()
    np.savetxt(ff2, data, fmt="%12.6G")
    fa = open(ff2.name)
    ff2.close()
    this = fa.readlines()
    f.writelines(this)
    f.close()

  def writeConf(self,confFile=None):
    '''
    write a configuration file

    Options:
      confFile : The name of the configuration file
                 If set to None, we use self.confFile
                 if thats not set we use 'config_files/sentinel0.conf'
  
    '''
    if not hasattr(self,'confTxt'):
      self.generateConfTxt()
    if not hasattr(self,'confFile'):
      self.confFile = 'config_files/sentinel0.conf'
    if self.confFile == None:
      self.confFile = 'config_files/sentinel0.conf'
    if confFile != None:
      self.confFile = confFile

    self.mkdir(self.confFile)
    if open(self.confFile,'w').write(self.confTxt) != None:
      self.error('error writing self.confTxt to configuration file %s'%self.confFile,fatal=True)

  def error(self,msg,fatal=False):
    '''
    Error reporting and exit if fatal set
    '''
    print >> sys.stderr, msg
    if fatal:
      sys.exit(1)

  def mkdir(self,name):
    '''
    Make directory
    '''
    try:
      os.makedirs(os.path.dirname(name))
    except OSError as exc: # Python >2.5
      if exc.errno == errno.EEXIST:
        pass
    else: raise

  def crossPlot(self,truth,est,params,filename=None):
    '''
    Cross plot

    given data from 'truth' and estimated ('est')
    plot the list of parameters params 
    for both with uncertainty (1.96 sd) 
    for the estimated

    The datasets can be filenames of PARAMETERS format files
    or dictionaries (e.g. those resulting from self.loadData(file)
    
    Options:
      filename  : output graph to file (e.g. plots/xx.png)
    '''
    if type(truth) == str:
      truth = self.loadData(truth)
    if type(est) == str:
       est = self.loadData(est)

    plt.ion()
    plt.clf()
    dd = np.in1d(truth['time'],est['time'])
    plt.plot([0.,1.],[0.,1.],'k-',label='1:1')
    for k in params:
      try:
        plt.errorbar(truth[k][dd],est[k],est['sd-' + k]*1.96)
        plt.plot(truth[k][dd],est[k],'*',label=k)
      except:
        pass
    plt.legend(loc='best')
    plt.xlabel('truth')
    plt.ylabel('est')
    if filename:
      self.mkdir(filename)
      plt.savefig(filename)
    else:
      plt.show()
  
  def paramPlot(self,truth,est,params,filename=None):
    '''
    Parameter plot

    given data from 'truth' and estimated ('est')
    plot the list of parameters params 
    for both, with uncertainty (1.96 sd) 
    for the estimated
    
    Options:
      filename  : output graph to file (e.g. plots/xx.png)
    '''
    plt.ion()
    plt.clf()
    npa = len(params)
    sq1 = int(np.sqrt(npa))
    sq2 = int(npa/float(sq1))
    fig = plt.figure(1)
    fig.clf()
    for jj in xrange(npa):
      k = params[jj]
      ax = fig.add_subplot(sq1,sq2,jj+1)
      ax.plot(truth['time'],truth[k])
      ax.errorbar(est['time'],est[k],est['sd-' + k]*1.96)
      ax.plot(est['time'],est[k],'ro',label=k)
      ax.legend(loc='best')
    if filename:
      self.mkdir(filename)
      plt.savefig(filename)
    else:
      plt.show()


  def generateConfTxt(self):
    '''
    Generate the default configuration file
    '''
    self.confTxt = '''
# configuration file
# for sentinel 


[parameter]
location = ['time']
limits = [[1,365,1]]
names=gamma_time,xlai, xhc,  rpl,  xkab, scenesc, xkw, xkm,   xleafn, xs1,xs2,xs3,xs4,lad
solve = [1]*len($parameter.names)
help_solve='flags for which state vector elements to solve for'

[parameter.result]
filename = 'output/rse1/rse1_test.dat'
help_filename="state vector results file"
format = 'PARAMETERS'

[parameter.assoc_solve]
gamma_time = 0
lad = 0
xs3 = 0
xs4 = 0
xhc = 0
rpl = 0
scenesc =  0
xs2 = 0

[parameter.x]
datatype = x
names = $parameter.names
default = [100,0.995,5,0.01,0.995,0.0,0.995,0.995,1.5,1.0,0,0,0,5]
help_default = "Set the parameter default values"
apply_grid = True
sd = [1.]*len($parameter.names)
bounds = [[0.01,0.99]]*len($parameter.names)
#state = data/rse1_init.dat
invtransform=$parameter.names
transform=$parameter.names

[parameter.x.assoc_transform]
xlai=np.exp(-xlai/2.)
xkab=np.exp(-xkab/100.)
#xkar=np.exp(-xkar/100.)
xkw=np.exp(-xkw*50.)
xkm=np.exp(-100.*xkm)

[parameter.x.assoc_invtransform]
xlai=-2.*np.log(xlai)
xkab=-100.*np.log(xkab)
#xkar=-100.*np.log(xkar)
xkw=-(1./50.)*np.log(xkw)
xkm=-(1./100.)*np.log(xkm)

[parameter.x.assoc_bounds]
gamma_time = 0.000001,100000
xlai = 0.067,0.995
xhc = 0.01,5.0
rpl = 0.001,0.10
xkab = 0.135,1
scenesc = 0.0,1
xkw = 0.135,1
xkm = 0.135,1
xleafn = 0.8,2.5
xs1 = 0.00, 1.5
xs2 = -2, 2
xs3 = -0.05,0.05
xs4 = -0.03,0.03
lad = 1,5 


[general]
is_spectral = True
calc_posterior_unc=True
help_calc_posterior_unc="Posterior uncertainty calculations"
write_results=True
doplot=True
help_doplot='plotting'
plotmod=30
help_plotmod='frequency of plotting'
plotmovie=False
epsilon=10e-15
help_epsilon="Epsilon"

[general.optimisation]
randomise=False

[operator]
prior.name=Operator
prior.datatypes = x,y
obs.name=Observation_Operator
obs.datatypes = x,y

[operator.prior.x]
names = $parameter.names[1:]
datatype = x

[operator.prior.y]
control = 'mask'.split()
names = $parameter.names[1:]
sd = [10000.0]*len($operator.prior.y.names)
help_sd='set the prior sd'
datatype = y
state = $parameter.x.default[1:]
help_state = "Set the prior state vector"

[operator.prior.y.result]
filename='output/rse1/rse1_test_prior.dat'
help_filename = 'prior filename'

[operator.obs.rt_model]
model=semidiscrete1
use_median=True
help_use_median = "Flag to state whether full bandpass function should be used or not. If True, then the median wavelength of the bandpass function is used"
bounds = [400,2500,1]
help_bounds = "The spectral bounds (min,max,step) for the operator'
ignore_derivative=False
help_ignore_derivative = "Set to True to override loading any defined derivative functions in the library and use numerical approximations instead"

[operator.obs.x]
names = $parameter.names[1:]
sd = [1.0]*len($operator.obs.x.names)
datatype = x

[operator.obs.y]
control = 'mask vza vaa sza saa'.split()
names = ['433-453','457.5-522.5','542.5-577.5','650-680','697.5-712.5','732.5-747.5','773-793','784.5-899.5','855-875','935-955','1565-1655','2100-2280']
sd = ["0.004", "0.00416142",  "0.00440183", "0.00476245", "0.00489983", "0.00502003","0.00516772", "0.00537035", "0.00544934", "0.0057241", "0.00800801","0.01" ]
datatype = y
state = 'data/rse1_test.100.dat'
help_state='set the obs state file'

[operator.obs.y.result]
filename = 'output/rse1/rse1_test_fwd.dat'
help_filename = 'forward modelling results file'
format = 'PARAMETERS'
''' 


def main():
  params = ['xlai','xkab','xkw','xkm','xleafn','xs1']

  s = Sentinel()
  # generate parameter dataset
  s.parameters(np.arange(1,366),'input/truth.dat')
  # generate noise free observations (every 5 days)
  s.fwdModel('input/truth.dat','input/sentinelClean.dat',every=5)
  s.fwdModel('input/truth.dat','input/sentinelGapClean.dat',every=5,prop=0.5,WINDOW=5)
  # add noise to the observations
  s.addNoiseToObservations('input/sentinelClean.dat','input/sentinel.dat',nMag=1.0)
  s.addNoiseToObservations('input/sentinelGapClean.dat','input/sentinelGap.dat',nMag=1.0)

  # solve for parameter estimate, using one date at a time
  # e.g. test for clean data starting from the correct value
  # This should work perfectly & go stright to the the solution
  # So, this is a sanity check for the solver mainly, but it is also
  # interesting to look at the distribution of uncertainties
  # in the plot plots/sentinelClean1_pplot.png
  # When we add noise, we expect the solution to lie somewhere in these bounds
  # generally.
  s.solveSingle('input/sentinelClean.dat','output/sentinelClean1.dat',initial='input/truth.dat')
  s.crossPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinelClean1.dat'),\
            params,filename='plots/sentinelClean1_xplot.png')
  s.paramPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinelClean1.dat'),\
            params,filename='plots/sentinelClean1_pplot.png')

  # e.g. test for clean data not starting from the correct value
  # Ideally, this would provide the same result, but that is unlikely
  # as the cost function is quite flat around the minimum.
  # But actually it does a pretty good job other than perhaps the 
  # first sample, which suggests that we need to reconsider
  # starting positions if we can afford it (e.g. go through the
  # series backwards as well) 
  s.solveSingle('input/sentinelClean.dat','output/sentinelClean2.dat')
  s.crossPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinelClean2.dat'),\
            params,filename='plots/sentinelClean2_xplot.png')
  s.paramPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinelClean2.dat'),\
            params,filename='plots/sentinelClean2_pplot.png')

  # e.g. test normal, i.e. data with noise. This is now a realistic test
  # of the solver for noisy data. We expect the result to be similar
  # to the truth (mostly within the 90% CI) 
  # Again, this should do a pretty decent job
  s.solveSingle('input/sentinel.dat','output/sentinel.dat')
  s.crossPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinel.dat'),\
            params,filename='plots/sentinel_xplot.png')
  s.paramPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinel.dat'),\
            params,filename='plots/sentinel_pplot.png')

  s.solveSingle('input/sentinelGap.dat','output/sentinelGap.dat')
  s.crossPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinelGap.dat'),\
            params,filename='plots/sentinelGap_xplot.png')
  s.paramPlot(s.loadData('input/truth.dat'),s.loadData('output/sentinelGap.dat'),\
            params,filename='plots/sentinelGap_pplot.png')

  # now, taking the result of the single inversion as the inital estimate (at the data points)
  # try to solve for all time samples using regularisation
  # If this takes some time, watch the world go by in the grpahs such as
  # in output/sentinel_O_1_gamma_100_resul*.png
  s.solveRegular('input/sentinel.dat','output/sentinel_O_1_gamma_150.dat',\
            modelOrder=1,gamma=150,initial='output/sentinel.dat')
  s.solveRegular('input/sentinelGap.dat','output/sentinelGap_O_1_gamma_150.dat',\
            modelOrder=1,gamma=150,initial='output/sentinelGap.dat')
  s.solveRegular('input/sentinel.dat','output/sentinel_O_2_gamma_530.dat',\
            modelOrder=2,gamma=530,initial='output/sentinel.dat')
  s.solveRegular('input/sentinelGap.dat','output/sentinelGap_O_2_gamma_530.dat',\
            modelOrder=2,gamma=530,initial='output/sentinelGap.dat')



if __name__ == "__main__":
  import sys
  sys.exit(main())

