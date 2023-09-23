from xspec import*
import matplotlib.pyplot as plt
from astropy.io import fits

list_gamma1 = []
list_gamma2 = []
list_A = []
list_E = []

list_time = []
list_time_err = []

list_A_dn = []
list_A_up = []

list_E_dn = []
list_E_up = []

list_gamma1_dn = []
list_gamma1_up = []

list_gamma2_dn = []
list_gamma2_up = []

a = 1

f_out = open('out2.txt', 'w')

line = ['n', 't1', 't2', 'gamma1', 'g1_', 'g1+', 'E', 'E_', 'E+', 'gamma2', 'g2_', 'g2+', 'A', 'A_', 'A+', 'chi2', 'dof']
tmp = line
f_out.write("{0:<6}{1:<15}{2:<15}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}{11:<10}{12:<10}{13:<10}"
            "{14:<10}{15:<20}{16:<20}\n".format(*tmp))


def build_model(model_name):
    #

    global a
    Fit.perform()
    Fit.nIterations = 500

    Fit.query = 'yes'
    Plot.device = 'xv'
    Plot.xAxis = "keV"
    # Plot("ldata delchi")
    Plot("uf delchi")

    xvals = Plot.x()
    y_model = Plot.model()
    xerros = Plot.xErr()
    yvals = Plot.y()
    yerros = Plot.yErr()

    b = len(y_model)

    y_err = []
    y_err0 = []

    for i in range(b):
        y_err0 += [0]

    for i in range(len(y_model)):
        y_err += [(yvals[i] - y_model[i])/yerros[i]]

    for i in range(b-len(y_model)):
        y_err += [0]
    for i in range(b-len(y_model)):
        y_model += [0]
    for i in range(b-len(y_model)):
        yvals += [0]

    plt.figure(figsize=(10, 5))
    plt.subplot(211)
    plt.title(spectrum_name)
    plt.step(xvals[:b], y_model[:b], label=r'model', c='r')
    plt.errorbar(xvals[:b], yvals[:b], xerr=xerros, yerr=yerros, fmt='None', c='0.4')
    plt.semilogx()
    plt.semilogy()
    plt.ylabel(r'Photons, ' + r'$s^{-1}cm^{-2}keV^{-1}$', fontsize=10)
    plt.xlabel(None)
    plt.legend(['model', 'data'], loc='best', fontsize=9)

    plt.subplot(212)
    plt.errorbar(xvals[:b], y_err[:b], xerr=xerros, yerr=0.7, fmt='None', c='0.4')
    plt.plot(xvals[:b], y_err0[:b], c='#00FF00')
    plt.semilogx()
    plt.ylabel(r'$\chi$', fontsize=10)
    plt.xlabel(r'Energy, $keV$', fontsize=10)

    plt.savefig('1_6spec_' + model_name + '_' + str(spectrum_name) + '.png')
    plt.show()


for i in range(60):

    spectrum_name = 'KW20030426_T03301_1_sp' + str(i+5) + '_gr10.pha'
    AllData.clear()

    s1 = Spectrum(spectrum_name)

    howmuch = Spectrum.dataGroup
    AllModels.lmod("grb")

    model_name1 = "npow"
    m1 = Model(model_name1)
    # m1.setPars(2.8, 2)
    build_model(model_name1)

    chi1 = Fit.statistic
    bins1 = Fit.dof + Fit.nVarPars

    dof1 = Fit.dof
    model_name2 = "nbknpow"
    m2 = Model(model_name2)
    build_model(model_name2)
    chi2 = Fit.statistic
    bins2 = Fit.dof + Fit.nVarPars

    dof2 = Fit.dof

    power1 = m2.nbknpow.PhoIndx1
    power_val1 = power1.values[0]

    Fit.error("maximum 100.")
    Fit.error("1. 1")
    power1_dn = power1.error[0] - power1.values[0]
    power1_up = power1.error[1] - power1.values[0]

    power1_dn = round(power1_dn, 3)
    power1_up = round(power1_up, 3)

    Energy = m2.nbknpow.BreakE

    Fit.error("1. 2")
    E_dn = Energy.error[0] - Energy.values[0]
    E_up = Energy.error[1] - Energy.values[0]

    E_dn = round(E_dn, 3)
    E_up = round(E_up, 3)

    power2 = m2.nbknpow.PhoIndx2
    power_val2 = power2.values[0]

    Fit.error("1. 4")
    power2_dn = power2.error[0] - power2.values[0]
    power2_up = power2.error[1] - power2.values[0]

    power2_dn = round(power2_dn, 3)
    power2_up = round(power2_up, 3)

    A = m2.nbknpow.norm
    A_val = A.values[0]

    Fit.error("1. 3")

    A_dn = A.error[0] - A.values[0]
    A_up = A.error[1] - A.values[0]

    A_dn = round(A_dn, 3)
    A_up = round(A_up, 3)

    statistic_value = Fit.ftest(chi2, dof2, chi1, dof1)

    if statistic_value < 0.001:
        print('nbknpow')
        gamma2 = power_val2
        if power2_dn != '-':
            gamma2_dn = power2_dn
            gamma2_up = power2_up
        gamma2 = round(gamma2, 3)

    else:
        print('npow')
        gamma2 = '--'
        energy = '--'
        gamma2_dn = '--'
        gamma2_up = '--'

    # create table
    hdul = fits.open(spectrum_name)

    line = [i + 5, hdul[1].header['TIME_OBS'], hdul[1].header['TIME_END'], round(power_val1, 3), power1_dn, power1_up, round(Energy.values[0], 3),
            E_dn, E_up, gamma2, power2_dn, power2_up, round(A_val, 3), A_dn, A_up, round(chi2, 3), dof2]
    list_gamma1 += [power_val1]

    if gamma2 != '--':
        list_gamma2 += [gamma2]

        list_gamma2_dn += [power2_dn]
        list_gamma2_up += [power2_up]
    else:
        list_gamma2 += [0]

        list_gamma2_dn += [0]
        list_gamma2_up += [0]

    list_A += [A_val]
    list_A_dn += [A_dn]
    list_A_up += [A_up]

    list_gamma1_dn = [power1_dn]
    list_gamma1_up = [power1_up]

    if Energy != '--':
        list_E += [round(Energy.values[0], 3)]
        list_E_dn += [E_dn]
        list_E_up += [E_up]
    else:
        list_E += [0]
        list_E_dn += [0]
        list_E_up += [0]
    list_time += [hdul[1].header['TSTART'] + hdul[1].header['T0TIME']]

    if len(list_time) >= 2:
        list_time_err += [(list_time[-1] - list_time[-2])/2]
    else:
        list_time_err += [0]

    tmp = line
    f_out.write("{0:<6}{1:<15}{2:<15}{3:<10}{4:<10}{5:<10}{6:<10}{7:<10}{8:<10}{9:<10}{10:<10}{11:<10}{12:<10}{13:<10}"
                "{14:<10}{15:<20}{16:<20}\n".format(*tmp))

    if i == 0:
        break
    # print(list_time)


plt.figure(figsize=(10, 5))

plt.subplot()
plt.errorbar(list_time, list_E, xerr=0.6, yerr=[list_E_dn, list_E_up], c='0.5', zorder=3, capsize=0, linestyle="None")

plt.ylabel(r'E', fontsize=10)
plt.xlabel(r'Time', fontsize=10)

plt.savefig('grafik2.png')
plt.show()
