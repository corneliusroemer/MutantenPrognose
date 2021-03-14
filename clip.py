import pandas as pd
from scipy.optimize import curve_fit
from scipy import stats
import datetime as dt
import numpy as np

df = pd.read_clipboard(parse_dates=True, index_col=0, dtype=int)
df
start_date = df.index[0]
(df.index - start_date).days
dfm = df
dfm = dfm.assign(t=lambda x: (
    x.index - start_date).days).assign(tlog=lambda x: np.log(x.t))
dfm = dfm.assign(nlog=lambda x: np.log(x['Impfung pro Tag']))[8:]
dfm
dfm.plot(x='tlog', y='nlog')

dfm = dfm[10:]
res = stats.linregress(dfm.tlog, dfm.nlog)
res.stderr

dfm = dfm.drop(columns='fit')
dfm = dfm.assign(fit=lambda x: np.exp(x.tlog*res.slope+res.intercept))
plot = dfm.plot(y=['Impfung pro Tag', 'fit'],
                xlim=(dfm.index[0], dfm.index[80]), figsize=(6.75, 3.5),
                title="Deutsche Impfrate im 7T-Mittel und Fit mit Potenzkurve",
                legend=False)

plot.set(ylabel="Tägliche Zahl an verimpften Dosen")
plot.get_figure().text(0.92, -0.08, "Datenstand: " + "21. Feb 2021, " +
                       "Datenquelle: RKI, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
fig = plot.get_figure().savefig("fit.png")

fig = plt.figure(num=None, figsize=(6.75, 3.5), facecolor='w', edgecolor='k')
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
ax.set_xlim(left=d.index[ps-1], right=d.index[-1]+pd.Timedelta(days=1))
ax.yaxis.set_major_locator(ticker.MultipleLocator(25))

plt.title("Inzidenzentwicklung & Prognose für Wildtyp & Mutanten in BaWü", size=11)
plt.ylabel("7-Tage Inzidenz")

plt.grid(axis='y', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.26, zorder=+1)

plt.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.15)
fig.text(0.94, 0.03, "Datenstand: " + "21. Feb 2021, " +
         "Datenquelle: Landesgesundheitsamt BaWü, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
plt.legend(prop={'size': 7})
plt.savefig('bawu.png', dpi=400)


ix = pd.date_range(start=dfm.index[0], end=dfm.index[-1]+dt.timedelta(100))
dfm = dfm.reindex(ix)
