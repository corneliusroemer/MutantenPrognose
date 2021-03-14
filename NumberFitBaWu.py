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


df = df.assign(week_cases=week_diff('Cases'))
df5=df
# %%
start = 8
end = 60
start_date = df5.index[start]
lvi = df5.index.get_loc(df5.week_cases.last_valid_index())
wt = df5[start:lvi].week_cases
x = range(len(wt))
print(wt.index[0])

px = np.linspace(-start, end, end+start)
# add trends like: lockerung from now onwards
# add seasonality (modulate b) with sine
# add vaccination (modulate b) with 14d delayed curve
amp =np.log(10)
ph = 30
lock = 0.02
lock_start = 46
mut_effect = 0.1  
vacc_const =  2.6/10/83/30*0.7
lock=0
amp=0
lock=0

def func(x, a, b, c):
    lock_effect = lock*np.heaviside(x-lock_start, 0.5)*(x-lock_start) + amp*(np.cos((x-ph)*(1/365)*math.pi*2)-1)
    lock_effect += np.exp(-vacc_const*(x**2))
    return a*np.exp((b)*x+lock_effect)+c*np.exp((b+mut_effect)*x+lock_effect)


def unfunc(x, a, b, c):
    lock_effect = lock*np.heaviside(x-lock_start, 0.5)*(x-lock_start) + amp*(unp.cos((x-ph)*(1/365)*math.pi*2)-1)
    lock_effect += np.exp(-vacc_const*(x**2))
    return a*unp.exp((b)*x+lock_effect)+c*unp.exp((b+mut_effect)*x+lock_effect)


popt_wt, pcov_wt = curve_fit(func, np.array(x), wt, p0=[140000, -0.03, 1400])

a, b, c = unc.correlated_values(popt_wt, pcov_wt)
b
py = (unfunc(px, a, b, c))

nom_wt2 = unp.nominal_values(py)
std_wt2 = unp.std_devs(py)

d = pd.DataFrame({'tot_fit': nom_wt2, 'tot_upper': nom_wt2+2*np.sqrt(std_wt2**2), 'tot_lower':
                  nom_wt2-2*np.sqrt(std_wt2**2)}, index=pd.date_range(start=df5.index[0], periods=end+start))
d = pd.concat([df5, d], axis=1)
d = d/110

d.describe()
#fit to current case numbers to make continuous
lvi = d.index.get_loc(d.week_cases.last_valid_index())
# d.tot_fit = d.tot_fit*d.week_cases[lvi]/d.tot_fit[lvi]
# d.to_csv('base+vacc+season.csv',columns=['week_cases','tot_fit'])

print(a, b, c)

ps = 1
lvi = d.index.get_loc(d.week_cases.last_valid_index())
fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
locale.setlocale(locale.LC_ALL, 'de_DE')
locale.setlocale(locale.LC_NUMERIC, 'de_DE')
plt.rcParams['axes.formatter.use_locale'] = True
ax.bar(d.index[ps:], d['week_cases'][ps:], color='#d94f45', label='Wildtyp')
ax.plot(d.index[ps:], d['tot_fit'][ps:], color='black',
        label='Inzidenz-Prognose mit 95% Konfidenzintervall')
ax.fill_between(d.index[lvi:], d['tot_lower'][lvi:],
                d['tot_upper'][lvi:], color='black', alpha=0.1, linewidth=0)


locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# locator = mdates.WeekdayLocator()
formatter = mdates.ConciseDateFormatter(locator, formats=[
                                        '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.set_xlim(left=d.index[ps-1], right=d.index[-1]+pd.Timedelta(days=1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(25))

plt.title(
    "Inzidenzentwicklung & Prognose für Wildtyp & Mutanten in Deutschland", size=11)
plt.ylabel("7-Tage Inzidenz")

plt.grid(axis='y', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.26, zorder=+1)

plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
fig.text(0.94, 0.03, "Datenstand: " + "24. Feb 2021, " +
         "Datenquelle: RKI & Risklayer , Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
plt.legend(prop={'size': 7})
plt.savefig('de.png', dpi=400)
# %%

# %%
# %%

# %%
