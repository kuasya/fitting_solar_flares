from xspec import*
import matplotlib.pyplot as plt

# print('1st task')

# Spectrum("file1.pha")
# Model("powerlaw")
# Fit.perform()
# Plot.device = "/xs"
# Plot("data")

s1 = Spectrum("KW20150507_T45695_1_sp1_4_gr10.pha")
# s1.show()
# Model.showList()
m1 = Model("powerlaw")
m1.setPars(2.8, 2)

Fit.perform()
Fit.nIterations = 500
print()
print('result')
# print(s1)

# Plot.device = "xserver-xorg"
# plot deviceplot type
# Plot("model")
# Plot.hardcopy('fig3.png') #no this attribute

# Plot.x()
# Plot.y()
# Plot.grid(True)
# Plot.legend(loc='best', fontsize=12)
# Plot.savefig('figure_with_legend.png')
# # Plot.yLog
# Plot.show()

Plot.setID('energy')
Plot.setGroup("1-**")
Plot.iplot("data delchi")
# log y
# lwidth 3
# time off
# hard fig2_py.png/png
# help(Plot)


