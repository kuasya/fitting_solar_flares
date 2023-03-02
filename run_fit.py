"""
Example of doing BXA in X-spec
"""
#import bxa.xspec as bxa
from xspec import *
from scipy.stats import chi2
import math

Fit.statMethod = 'chi2'
Fit.nIterations = 500
Plot.xAxis = 'keV'

AllModels.lmod('grb')
m=Model("npow+gaussian+nuclear+gaussian")

# set model parameters ranges
#                         val, delta, min, bottom, top, max
#m.npow.PhoIndex.values = "3,0.01,-3,-2,10,10" 
m.npow.PhoIndex.values = "3,0.01,0,0,10,10" 
m.npow.enorm.values = "500,-1,1,1,1e3,1e5" 
#m.npow.norm.values = "1,0.01,1e-10,1e-10,1e+2,1e+2" 
m.npow.norm.values = "1,0.01,0,0,100,100" 

m.gaussian.LineE.values = "511,-1,1,1,1e+20,1e+24"
m.gaussian.Sigma.values = "3,-1,0.1,0.1,100,200"
#m.gaussian.norm.values = "1,0.01,1e-10,1e-10,1e+2,1e+2"
m.gaussian.norm.values = "1,0.01,0,0,100,100"

#m.nuclear.norm.values = "1,0.01,1e-10,1e-10,1e+2,1e+2"
m.nuclear.norm.values = "1,0.01,0,0,100,100"
#m.pion50.norm.values = "1,0.01,1e-10,1e-10,1e+20,1e+24"

m.gaussian_4.LineE.values = "2223,-1,1,1,1e+20,1e+24"
m.gaussian_4.Sigma.values = "0.1,-1,0.1,0.1,100,200"
#m.gaussian.norm.values = "1,0.01,1e-10,1e-10,1e+2,1e+2"
m.gaussian_4.norm.values = "1,0.01,0,0,100,100"



flare='KW20170906_T42923'
s_range=['1','2']

spectra=['26_64','26_56', '57','58', '59','60', '61','62', '63','64']
#spectra=['60']

out_file='KW20170906_T42923_npow_511gauss_nucl_22gauss.dat'
f=open(out_file,'w')

str1="PHAname\tpow\tpow_dn\tpow_up\tA\tA_dn\tA_up\t511flux\t511flux_dn\t511flux_up\tnuclear\tnuclear_dn\tnuclear_up\t2.2flux\t2.2flux_dn\t2.2flux_up\tchi2\tdof\tprob\n"
f.write(str1)
keV511_uplims=[]
MeV22_uplims=[]

for spec in (spectra):
    
    AllData.clear()
    spec_name_1=flare+'_'+s_range[0]+'_sp'+spec+'.pha'
    print "Analyzing spectrum", spec_name_1
    AllData+=spec_name_1
    spec_name_2=flare+'_'+s_range[1]+'_sp'+spec+'.pha'
    AllData+=spec_name_2
    
    Xset.logChatter = 10
    log_name=flare+'_'+'_sp'+spec+'.log'
    logFile = Xset.openLog(log_name)


    Fit.nIterations = 500
    Fit.query = 'yes'
        
    AllData.ignore("bad")
    AllData.ignore("1:**-200. 420.-**")
    
    Fit.perform()
    chi2_=Fit.statistic
    dof=Fit.dof
        
    power=m.npow.PhoIndex
    power_val=power.values[0]
    Fit.error("2.7 1")
    print 'power=',power.values
    print 'power.error[0]=',power.error[0]
    power_dn=power.error[0]-power.values[0]
    power_up=power.error[1]-power.values[0]
    """
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "PhoIndex statistics"
    print power.values[0],power.error[0],power.error[1]
    lst_stat=ChainManager.stat(1)
    print lst_stat
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    """
    A=m.npow.norm
    A_val=A.values[0]
    Fit.error("2.7 3")
    A_dn=A.error[0]-A.values[0]
    A_up=A.error[1]-A.values[0]
    """
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "A statistics"
    print A.values[0],A.error[0],A.error[1]
    print ChainManager.stat(c,3)    
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    """
    gauss511=m.gaussian.norm
    gauss511_val=gauss511.values[0]
    Fit.error("2.7 6")
    gauss511_dn=gauss511.error[0]-gauss511.values[0]
    gauss511_up=gauss511.error[1]-gauss511.values[0]
    keV511_uplims.append(gauss511.error[1])
    
    nuclear=m.nuclear.norm
    nuclear_val=nuclear.values[0]
    Fit.error("2.7 7")
    nuclear_dn=nuclear.error[0]-nuclear.values[0]
    nuclear_up=nuclear.error[1]-nuclear.values[0]
    """
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "Nuckear Flux statistics"
    print nuclear.values[0],nuclear.error[0],nuclear.error[1]
    print ChainManager.stat(4)    
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    """    
    gauss22=m.gaussian_4.norm
    gauss22_val=gauss22.values[0]
    Fit.error("2.7 10")
    gauss22_dn=gauss22.error[0]-gauss22.values[0]
    gauss22_up=gauss22.error[1]-gauss22.values[0]
    MeV22_uplims.append(gauss22.error[1])
    """
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "2.2 MeV line statistics"
    print gauss22.values[0],gauss22.error[0],gauss22.error[1]
    print ChainManager.stat(7)    
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    """    
    prob = 1.0 - chi2.cdf(Fit.statistic, Fit.dof)    
    out_str="%s\t%.3f\t%.3f\t%.3f\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.2f\t%d\t%.2f\n" % \
        (spec,power_val,power_dn,power_up,A_val,A_dn,A_up,gauss511_val,gauss511_dn,gauss511_up, \
            nuclear_val,nuclear_dn,nuclear_up,gauss22_val,gauss22_dn,gauss22_up,chi2_,dof,prob)
        
    f.write(out_str)
    print keV511_uplims
    print MeV22_uplims
    
f.close()
    

