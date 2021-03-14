# %%

import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import uncertainties as unc
import uncertainties.unumpy as unp
import matplotlib.dates as mdates
import locale
import datetime as dt
import matplotlib.ticker as ticker
import dateutil.parser
import math
# %%

df = pd.read_csv('bawu.csv', index_col=0, parse_dates=True)

df['Cases'].rolling(window=8).apply(lambda x: x[-1]-x[0])


def week_diff(column):
    return lambda x: x[column].rolling(window=8).apply(lambda x: x[-1]-x[0])


df = df.assign(week_cases_bw=week_diff('Cases'))

# %%
df3 = pd.read_csv('risklayerde.csv', index_col=0,
                  parse_dates=True, dayfirst=True)

df3['cases'].rolling(window=8).apply(lambda x: x[-1]-x[0])


def week_diff(column):
    return lambda x: x[column].rolling(window=8).apply(lambda x: x[-1]-x[0])


# %%
df3 = df3.assign(week_cases_de=week_diff('cases'))

# %%
df4 = df3.join(df).tail(25)
# %%
df5 = df4.assign(inz_de_ex_bw=lambda x: (
    x.week_cases_de - x.week_cases_bw)/(830-110.7))
df5 = df5.assign(inz_bw=lambda x: x.week_cases_bw/110.7)
# %%
df5
# %%
df5.plot(y=['inz_de_ex_bw', 'inz_bw'])
# %%
# Why not add some derived formula to calculate wildtype
# %%
fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
locale.setlocale(locale.LC_ALL, 'de_DE')
locale.setlocale(locale.LC_NUMERIC, 'de_DE')
plt.rcParams['axes.formatter.use_locale'] = True
ax.plot(df5.index, df5.inz_bw, color='b', label='Inzidenz Wildtyp')

locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# locator = mdates.WeekdayLocator()
formatter = mdates.ConciseDateFormatter(locator, formats=[
                                        '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)

plt.title("Inzidenz des Wildtyps in Baden-Württemberg", size=11)
plt.ylabel("7-Tage Inzidenz")
left, right = plt.xlim()
plt.xlim((left+1,right-1))
plt.grid(axis='y', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.26, zorder=+1)
plt.axvline(dt.datetime(2021,2,11),linestyle='dotted',label="Aufhebung Ausgangssperre")

plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
fig.text(0.94, 0.03, "Datenstand: " + "27. Feb 2021, " +
         "Datenquelle: Lagebericht LGA BaWü, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
plt.legend(prop={'size': 9})
plt.savefig('bawu_wildtyp.png', dpi=400)
# %%
df5=df.assign(inz_bw=lambda x: x.wt_fit/107)
# %%

df5=df.assign(R=lambda x:x.wt_fit.rolling(window=8).apply(lambda x: x[-1]/x[0]))
plt.plot(df5.R)
# %%
