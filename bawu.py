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



#%%
df = pd.read_csv('bawu.csv', index_col=0, parse_dates=True)

df['Cases'].rolling(window=8).apply(lambda x: x[-1]-x[0])


def week_diff(column):
    return lambda x: x[column].rolling(window=8).apply(lambda x: x[-1]-x[0])


df = df.assign(week_cases=week_diff('Cases'))
df = df.assign(week_en=week_diff('EN'))
df = df.assign(week_sa=week_diff('SA'))
df = df.assign(mutants=lambda x: x.week_en+x.week_sa)
df = df.assign(wt=lambda x: x.week_cases - x.mutants)

df[7:].plot(y="week_en")

start = 7
end = 56
start_date = df.index[start]
wt = df[start:].wt
mut = df[start:].mutants
x = range(len(wt))
wt.index[0]


def func(x, a, b):
    return a+b*x


popt_wt, pcov_wt = curve_fit(func, x, np.log(wt))
popt_mut, pcov_mut = curve_fit(func, x, np.log(mut))

a, b = unc.correlated_values(popt_wt, pcov_wt)
b
px = np.linspace(-start, end, end+start)
py = unp.exp(func(px, a, b))

nom_wt2 = unp.nominal_values(py)
std_wt2 = unp.std_devs(py)

a, b = unc.correlated_values(popt_mut, pcov_mut)

px = np.linspace(-start, end, end+start)
py = unp.exp(func(px, a, b))

nom_mut2 = unp.nominal_values(py)
std_mut2 = unp.std_devs(py)

d = pd.DataFrame({'nom': nom_mut2+nom_wt2, 'upper': nom_mut2+nom_wt2+2*np.sqrt(std_mut2**2+std_wt2**2), 'lower': nom_mut2 +
                  nom_wt2-2*np.sqrt(std_mut2**2+std_wt2**2)}, index=pd.date_range(start=df.index[0], periods=end+start))
d2 = pd.DataFrame({'mut': nom_mut2, 'mut_upper': nom_mut2+2*std_mut2, 'mut_lower': nom_mut2 -
                   2*std_mut2}, index=pd.date_range(start=df.index[0], periods=end+start))
d = pd.concat([df, d, d2], axis=1)
d = d/111

d.describe()

# %%
loi = len(df)
ps = 7
lvi = d.index.get_loc(d.week_cases.last_valid_index())
fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
locale.setlocale(locale.LC_ALL, 'de_DE')
locale.setlocale(locale.LC_NUMERIC, 'de_DE')
plt.rcParams['axes.formatter.use_locale'] = True
ax.bar(d.index[ps:], d['week_cases'][ps:], color='#d94f45', label='Wildtyp')
ax.bar(d.index[ps:], d['mutants'][ps:], color='#5285ed',
       label='Mutanten (B.1.1.7 & B.1.351)')
ax.plot(d.index[ps:], d['nom'][ps:], color='black',
        label='Inzidenz-Prognose mit 95% Konfidenzintervall')
ax.plot(d.index[lvi+1:], d['mut'][lvi+1:], color='#5285ed',
        label='Mutanten-Prognose mit 95% Konfidenzintervall')
ax.fill_between(d.index[lvi:], d['lower'][lvi:],
                d['upper'][lvi:], color='black', alpha=0.1, linewidth=0)
ax.fill_between(d.index[lvi:], d['mut_lower'][lvi:],
                d['mut_upper'][lvi:], color='#5285ed', alpha=0.1, linewidth=0)


locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# locator = mdates.WeekdayLocator()
formatter = mdates.ConciseDateFormatter(locator, formats=[
                                        '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.set_xlim(left=d.index[ps]) #right=d.index[-1]+pd.Timedelta(days=1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(25))

plt.title("Inzidenzentwicklung & Prognose für Wildtyp & Mutanten in BaWü", size=11)
plt.ylabel("7-Tage Inzidenz")

plt.grid(axis='y', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.26, zorder=+1)

plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
fig.text(0.94, 0.03, "Datenstand: " + "26. Feb 2021, " +
         "Datenquelle: Landesgesundheitsamt BaWü, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
plt.legend(prop={'size': 7},loc='upper right')
plt.savefig('bawu.png', dpi=400)

# %%
# %%
