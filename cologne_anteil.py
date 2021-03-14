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

df = pd.read_csv('cologne.csv', index_col=0, parse_dates=True)

df['Cases'].rolling(window=8).apply(lambda x: x[-1]-x[0])


def week_diff(column):
    return lambda x: x[column].rolling(window=8).apply(lambda x: x[-1]-x[0])


df = df.assign(week_cases=week_diff('Cases'))
df = df.assign(week_en=week_diff('EN'))
df = df.assign(week_sa=week_diff('SA'))
df = df.assign(mutants=lambda x: x.week_en+x.week_sa)
df = df.assign(wt=lambda x: x.week_cases - x.mutants)
df = df.assign(Anteil=lambda x: x.mutants/x.week_cases)
# %%
df = df.loc[df.index >= dateutil.parser.parse("2021-01-28")]
df

# %%
df = df.assign(ratio=lambda x: (x.Anteil/(1-x.Anteil)))
df = df.assign(log_ratio=lambda x: np.log(x.ratio))
df = df.assign(days=lambda x: (df.index - dt.datetime(2021, 1, 1)).days)
df
# %%
df.plot(y="log_ratio")
# %%
end = 100


def func(x, a, b):
    return a+b*x


popt, pcovt = curve_fit(func, df.days, df.log_ratio)

a, b = unc.correlated_values(popt, pcovt)
print(a, b)
px = np.linspace(1, end, end)
fit = unp.exp(func(px, a, b))
py = (fit/(fit+1))

nom = unp.nominal_values(py)
std = unp.std_devs(py)
# %%
nom
std
# %%

df2 = pd.DataFrame({'fit': nom, 'err': std}, index=pd.date_range(
    '2021-01-01', freq='D', periods=end))
# %%
plt.plot(px, nom)
plt.plot(px, nom+std)
plt.plot(px, nom-std)
plt.plot(df.days, df.Anteil, linestyle="None", marker='s')
# %%
df2.std
# %%
fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
locale.setlocale(locale.LC_ALL, 'de_DE')
locale.setlocale(locale.LC_NUMERIC, 'de_DE')
plt.rcParams['axes.formatter.use_locale'] = True
ax.plot(df2.index, df2.fit, color='black',
        label='VOC-Anteil-Prognose mit 95% Konfidenzintervall')
ax.fill_between(df2.index, df2.fit-2*df2.err,
                df2.fit+2*df2.err, alpha=0.2, linewidth=0)
ax.plot(df.index, df.Anteil, marker='+', markersize=10, color='r',
        linestyle="None", label="Gemessener VOC-Anteil (VOC=Varianten)")

locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# locator = mdates.WeekdayLocator()
formatter = mdates.ConciseDateFormatter(locator, formats=[
                                        '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.set_xlim(left=d.index[ps-1], right=d.index[-1]+pd.Timedelta(days=1))
ax.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
plt.title(
    "Entwicklung des Variantenanteils (B.1.1.7 & Co) in Köln", size=11)
plt.ylabel("Anteil der Variante an allen Fällen")

plt.grid(axis='y', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.26, zorder=+1)

plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
fig.text(0.94, 0.03, "Datenstand: " + "1.3.2021, " +
         "Datenquelle: Stadt @Koeln Twitter, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
plt.legend(prop={'size': 8})
plt.savefig('cologne_voc.png', dpi=400)

# %%
df2[50:70]
# %%
