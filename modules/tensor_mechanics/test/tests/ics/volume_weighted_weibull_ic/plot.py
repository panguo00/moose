#!/usr/bin/env python
import pandas
import matplotlib.pyplot as plt

ax = plt.gca()
histo50 = pandas.read_csv('histo50.csv')
histo50.plot(ax=ax,x='u_aux',y='n',label='50')

histo100 = pandas.read_csv('histo100.csv')
histo100.plot(ax=ax,x='u_aux',y='n',label='100')

histo500 = pandas.read_csv('histo500.csv')
histo500.plot(ax=ax,x='u_aux',y='n',label='500')

plt.savefig('weibull.pdf')
