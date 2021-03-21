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
df = pd.read_csv('vocde2.csv', index_col=0, parse_dates=True)
week_series = pd.date_range(start='2021-1-28', periods=len(df), freq='14D')
df = df.set_index(week_series)
# %%
# df = pd.read_csv('vocde.csv', index_col=0, parse_dates=True)
# week_series = pd.date_range(start='2021-1-14', periods=len(df), freq='7D')
# df = df.set_index(week_series)
# %%
df = df.assign(ratio=lambda x: (x.Anteil/(1-x.Anteil)))
df = df.assign(log_ratio=lambda x: np.log(x.ratio))
df = df.assign(days=lambda x: (df.index - dt.datetime(2020, 12, 31)).days)
# %%
df
# %%
df.plot(y="log_ratio")
# %%
end = 200


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
plt.plot(px, nom)
plt.plot(px, nom+std)
plt.plot(px, nom-std)
plt.plot(df.days, df.Anteil, linestyle="None", marker='s')

# %%
df2 = pd.DataFrame({'fit': nom, 'std': std}, index=pd.date_range(
    '2021-01-01', freq='D', periods=end))
df2 = df2.join(df)
df4 = df2
# %%
plt.plot(nom)

# # %%
# fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
# ax = fig.add_subplot(1, 1, 1)
# locale.setlocale(locale.LC_ALL, 'de_DE')
# locale.setlocale(locale.LC_NUMERIC, 'de_DE')
# plt.rcParams['axes.formatter.use_locale'] = True
# ax.plot(df4.index, df4.fit, color='black',
#         label='VOC-Anteil-Prognose mit 95% Konfidenzintervall')
# ax.fill_between(df4.index, df4.fit-2*df4['std'],
#                 df4.fit+2*df4['std'], alpha=0.2, linewidth=0)
# ax.plot(df.index, df.Anteil, marker='+', markersize=10, color='r',
#         linestyle="None", label="Gemessener VOC-Anteil (VOC=Varianten)")

# locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# # locator = mdates.WeekdayLocator()
# formatter = mdates.ConciseDateFormatter(locator, formats=[
#                                         '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
# ax.xaxis.set_major_locator(locator)
# ax.xaxis.set_major_formatter(formatter)
# ax.set_xlim(left=d.index[ps-1], right=d.index[-1]+pd.Timedelta(days=1))
# ax.yaxis.set_major_formatter(ticker.PercentFormatter(1.0))
# plt.title(
#     "Entwicklung des Variantenanteils (B.1.1.7 & Co) in Deutschland", size=11)
# plt.ylabel("Anteil der Variante an allen Fällen")

# plt.grid(axis='y', b=True, which='both', color='#999999',
#          linestyle='-', alpha=0.26, zorder=+1)

# plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
# fig.text(0.94, 0.03, "Datenstand: " + "3.3.21, " +
#          "Datenquelle: RKI-Varianten-Bericht 3.3.21, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
# plt.legend(prop={'size': 8})
# plt.savefig('vocde.png', dpi=400)
# # %%
len(nom)
# %%
df2 = pd.DataFrame({'fit': nom, 'std': std}, index=pd.date_range(
    '2021-01-01', freq='D', periods=end))
# %%
df4 = df.join(df2, how='right')
# %%
df4 = df.join(df2, how='right')
# %%

df3 = pd.read_csv('risklayerde.csv', index_col=0,
                  parse_dates=True, dayfirst=True)

df3['cases'].rolling(window=8).apply(lambda x: x[-1]-x[0])


def week_diff(column):
    return lambda x: x[column].rolling(window=8).apply(lambda x: x[-1]-x[0])


# %%
print(unp.exp(b), unp.exp(4*b), unp.exp(7*b))
# %%
df3 = df3.assign(week_cases=week_diff('cases'))
dfx = df3/830
dfx.to_csv('weekly_incidence.csv', columns=['week_cases'])
dfx.plot()
# %%
df3 = df3[df3.index >= dateutil.parser.parse("2021-01-01")]
df3
# %%
df3.week_cases.plot()

# %%
df5 = df4.join(df3, how='left')
df5 = df5.assign(mut=lambda x: x.fit*x.week_cases)
df5 = df5.assign(mut_std=lambda x: x['std']*x.week_cases)
df5 = df5.assign(wt=lambda x: (1-x.fit)*x.week_cases)
df5 = df5.assign(wt_std=lambda x: x['std']*x.week_cases)

# %%

ps = 15
lvi = df5.index.get_loc(df5['wt'].last_valid_index())
# fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
# ax = fig.add_subplot(1, 1, 1)
# locale.setlocale(locale.LC_ALL, 'de_DE')
# locale.setlocale(locale.LC_NUMERIC, 'de_DE')
# plt.rcParams['axes.formatter.use_locale'] = True
# ax.bar(df5.index, df5.wt/830, color='b', label='Wildtyp')
# ax.bar(df5.index, df5.week_cases/830, color='r', alpha=0.15, label='Varianten')
# locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# # locator = mdates.WeekdayLocator()
# formatter = mdates.ConciseDateFormatter(locator, formats=[
#                                         '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
# ax.xaxis.set_major_locator(locator)
# ax.xaxis.set_major_formatter(formatter)
# ax.set_xlim(left=df5.index[ps-1], right=d.index[lvi]+pd.Timedelta(days=1))
# plt.title(
#     "Entwicklung der Wildtyp-Inzidenz in Deutschland", size=11)
# plt.ylabel("7T Inzidenz pro 100k EW")

# plt.legend(prop={'size': 9})
# plt.grid(axis='y', b=True, which='both', color='#999999',
#          linestyle='-', alpha=0.26, zorder=+1)

# plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
# fig.text(0.94, 0.03, "Datenstand: " + "3.3.21, " +
#          "Datenquelle: Risklayer, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
# plt.savefig('wildtypde.png', dpi=400)
# %%
df5.wt.plot()
plt.plot(df5.mut[20:])
plt.plot(df5.wt[20:])
plt.yscale('log')
# %%

# ps = 15

# start = 10
# end = 100
# start_date = df5.index[start]
# lvi = df5.index.get_loc(df5.wt.last_valid_index())
# wt = df5[start:lvi].wt
# mut = df5[start:lvi].mut
# x = range(len(wt))
# wt.index[0]


# def func(x, a, b):
#     return a+b*x


# popt_wt, pcov_wt = curve_fit(func, x, np.log(wt))
# popt_mut, pcov_mut = curve_fit(func, x, np.log(mut))

# a, b = unc.correlated_values(popt_wt, pcov_wt)
# b
# px = np.linspace(-start, end, end+start)
# py = unp.exp(func(px, a, b))

# nom_wt2 = unp.nominal_values(py)
# std_wt2 = unp.std_devs(py)

# a, b = unc.correlated_values(popt_mut, pcov_mut)

# px = np.linspace(-start, end, end+start)
# py = unp.exp(func(px, a, b))

# nom_mut2 = unp.nominal_values(py)
# std_mut2 = unp.std_devs(py)

# d = pd.DataFrame({'tot_fit': nom_mut2+nom_wt2, 'tot_upper': nom_mut2+nom_wt2+2*np.sqrt(std_mut2**2+std_wt2**2), 'tot_lower': nom_mut2 +
#                   nom_wt2-2*np.sqrt(std_mut2**2+std_wt2**2)}, index=pd.date_range(start=df5.index[0], periods=end+start))
# d2 = pd.DataFrame({'mut_fit': nom_mut2, 'mut_upper': nom_mut2+2*std_mut2, 'mut_lower': nom_mut2 -
#                    2*std_mut2}, index=pd.date_range(start=df5.index[0], periods=end+start))
# d = pd.concat([df5, d, d2], axis=1)
# d = d/830

# # d.describe()
# # %%

# ps = 15
# lvi = d.index.get_loc(d.week_cases.last_valid_index())
# fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
# ax = fig.add_subplot(1, 1, 1)
# locale.setlocale(locale.LC_ALL, 'de_DE')
# locale.setlocale(locale.LC_NUMERIC, 'de_DE')
# plt.rcParams['axes.formatter.use_locale'] = True
# ax.bar(d.index[ps:], d['week_cases'][ps:], color='#d94f45', label='Wildtyp')
# ax.bar(d.index[ps:], d['mut'][ps:], color='#5285ed',
#        label='Mutanten (B.1.1.7 & B.1.351)')
# ax.plot(d.index[ps:], d['tot_fit'][ps:], color='black',
#         label='Inzidenz-Prognose mit 95% Konfidenzintervall')
# ax.plot(d.index[lvi+1:], d['mut_fit'][lvi+1:], color='#5285ed',
#         label='Mutanten-Prognose mit 95% Konfidenzintervall')
# ax.fill_between(d.index[lvi:], d['tot_lower'][lvi:],
#                 d['tot_upper'][lvi:], color='black', alpha=0.1, linewidth=0)
# ax.fill_between(d.index[lvi:], d['mut_lower'][lvi:],
#                 d['mut_upper'][lvi:], color='#5285ed', alpha=0.1, linewidth=0)


# locator = mdates.AutoDateLocator(minticks=6, maxticks=20)
# # locator = mdates.WeekdayLocator()
# formatter = mdates.ConciseDateFormatter(locator, formats=[
#                                         '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
# ax.xaxis.set_major_locator(locator)
# ax.xaxis.set_major_formatter(formatter)
# ax.set_xlim(left=d.index[ps-1], right=d.index[-1]+pd.Timedelta(days=1))
# ax.yaxis.set_major_locator(ticker.MultipleLocator(25))

# plt.title(
#     "Inzidenzentwicklung & Prognose für Wildtyp & Mutanten in Deutschland", size=11)
# plt.ylabel("7-Tage Inzidenz")

# plt.grid(axis='y', b=True, which='both', color='#999999',
#          linestyle='-', alpha=0.26, zorder=+1)

# plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
# fig.text(0.94, 0.03, "Datenstand: " + "24. Feb 2021, " +
#          "Datenquelle: RKI & Risklayer , Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
# plt.legend(prop={'size': 7})
# plt.savefig('de.png', dpi=400)
# %%
df5 = df5.assign(day=lambda x: (
    x.index-dateutil.parser.parse('20210101')).days)
# %%


class Scenario:
    pass


# Best Guess
bg = Scenario()
bg.amp = np.log(800)/2
bg.ph = 60
bg.lock = 0.025
bg.lock_start = 44
bg.mut_effect = 0.088
bg.vacc_const = 6/10/83/30*0.7
bg.lock2 = 0.015
bg.lock2_start = 68
bg.lock3 = -0.010
bg.lock3_start = 78

# Super Optimistic Seasonality and Vaccination
bot = Scenario()
bot.amp = np.log(3000)/2
bot.ph = 60
bot.lock = 0.025
bot.lock_start = 44
bot.mut_effect = 0.084
bot.vacc_const = 6/10/83/30*0.7
bot.lock2 = 0.015
bot.lock2_start = 68
bot.lock3 = -0.010
bot.lock3_start = 76

# Pessimistic
top = Scenario()
top.amp = np.log(300)/2
top.ph = 60
top.lock = 0.025
top.lock_start = 44
top.mut_effect = 0.092
top.vacc_const = 5/10/83/30*0.7
top.lock2 = 0.020
top.lock2_start = 68
top.lock3 = -0.014
top.lock3_start = 78

scenarios = {'best': bg, 'bottom': bot, 'top': top}
# %%
lvi = df5.index.get_loc(df5.wt.last_valid_index())
end = lvi
start = lvi-50
start_date = df5.index[start]
plot_end = lvi+28
wt = df5[start:lvi+1].week_cases
x = range(len(wt))
print(wt.index[0])

px = np.linspace(-start, end, end+start)
# add trends like: lockerung from now onwards
# add seasonality (modulate b) with sine
# add vaccination (modulate b) with 14d delayed curve

# %%

for k, v in scenarios.items():
    def func(x, a, b, c):
        lock_effect = v.lock3 * \
            np.heaviside(x-v.lock3_start, 0.5)*(x-v.lock3_start)
        lock_effect += v.lock2 * \
            np.heaviside(x-v.lock2_start, 0.5)*(x-v.lock2_start)
        lock_effect += v.lock*np.heaviside(x-v.lock_start, 0.5) * \
            (x-v.lock_start) + v.amp*(np.cos((x-v.ph)*(1/365)*math.pi*2)-1)
        lock_effect += np.exp(-v.vacc_const*(x**2))
        return a*np.exp((b)*x+lock_effect)+c*np.exp((b+v.mut_effect)*x+lock_effect)

    popt_wt, pcov_wt = curve_fit(
        func, df5.day[start:end+1], df5.week_cases[start:end+1], p0=[140000, -0.03, 1400])

    kwargs = {k: lambda x: func(x.day, *popt_wt)/830}
    df5 = df5.assign(**kwargs)
    print(popt_wt)

# %%
df5 = df5.assign(inz=df5.week_cases/830)
df5 = df5.assign(mutanten=df5.mut/830)
ratio = 0.5
best = df5.best[lvi:plot_end]
top = df5.top[lvi:plot_end]
bot = df5.bottom[lvi:plot_end]
df5 = df5.assign(t95=(top)*1.02)
df5 = df5.assign(b95=(bot)*0.98)
df5 = df5.assign(b50=((1-ratio)*best + ratio * bot)*0.98)
df5 = df5.assign(t50=((1-ratio)*best + ratio * top)*1.02)
t95 = df5.t95[lvi:plot_end]
b95 = df5.b95[lvi:plot_end]
t50 = df5.t50[lvi:plot_end]
b50 = df5.b50[lvi:plot_end]
df5.to_csv('scenarios.csv', columns=[
           'inz', 'mutanten', 'best', 't95', 'b95', 't50', 'b50'])
# %%

# fit to current case numbers to make continuous
# lvi = d.index.get_loc(d.week_cases.last_valid_index())
# d.tot_fit = d.tot_fit*d.week_cases[lvi]/d.tot_fit[lvi]
# d.to_csv('base+vacc+season.csv',columns=['week_cases','tot_fit'])

print(*popt_wt)

lvi = df5.index.get_loc(df5.wt.last_valid_index())
ps = lvi-50
# lvi = d.index.get_loc(d.week_cases.last_valid_index())
fig = plt.figure(num=None, figsize=(6.2, 3.5), facecolor='w', edgecolor='k')
ax = fig.add_subplot(1, 1, 1)
locale.setlocale(locale.LC_ALL, 'de_DE')
locale.setlocale(locale.LC_NUMERIC, 'de_DE')
plt.rcParams['axes.formatter.use_locale'] = True
# ax.fill_between(df5.index[lvi:plot_end], ((1-ratio)*best + ratio * bot-1)*0.98,
#                 ((1-ratio)*best + ratio * top+1)*1.02, color='black', alpha=0.2, linewidth=0, label='50% Konfidenzintervall')
# ax.fill_between(df5.index[lvi:plot_end], (bot-2)*0.97,
#                 (top+2)*1.03, color='black', alpha=0.1, linewidth=0, label='95% Konfidenzintervall')
ax.fill_between(df5.index[lvi:plot_end], b50,
                t50, color='black', alpha=0.2, linewidth=0, label='50% Konfidenzintervall')
ax.fill_between(df5.index[lvi:plot_end], b95,
                t95, color='black', alpha=0.1, linewidth=0, label='95% Konfidenzintervall')
ax.bar(df5.index[ps:plot_end], df5['week_cases']
       [ps:plot_end]/830, color='#d94f45', label='Wildtyp')
ax.bar(df5.index[ps:plot_end], df5['mut'][ps:plot_end]/830,
       color='#5285ed', label='Mutanten (B.1.1.7 & B.1.351)')
ax.plot(df5.index[ps:plot_end], df5.best[ps:plot_end], color='black', lw=0.5)
ax.plot(df5.index[lvi:plot_end], df5.best[lvi:plot_end],
        color='black', label='Inzidenz-Prognose')
print(df5.describe())

locator = mdates.AutoDateLocator(minticks=6, maxticks=15)
# locator = mdates.WeekdayLocator()
formatter = mdates.ConciseDateFormatter(locator, formats=[
    '%Y', '%-d.%b', '%-d', '%H:%M', '%H:%M', '%S.%f'], show_offset=False)
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.set_xlim(left=df5.index[ps], right=df5.index[plot_end])
# ax.yaxis.set_major_locator(ticker.MultipleLocator(25))

plt.title(
    "Inzidenzentwicklung & Prognose für Wildtyp & Mutanten in Deutschland", size=11, x=0.46)
plt.ylabel("7-Tage Inzidenz")

plt.grid(axis='y', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.26, zorder=-10)
plt.grid(axis='x', b=True, which='both', color='#999999',
         linestyle='-', alpha=0.06, zorder=-10)

plt.subplots_adjust(left=0.13, right=0.95, top=0.9, bottom=0.15)
fig.text(0.94, 0.03, "Datenstand: " + "21. März 2021, " +
         "Datenquelle: RKI & Risklayer, Analyse: @CorneliusRoemer", size=7, va="bottom", ha="right")
plt.legend(prop={'size': 8})
plt.savefig('de.png', dpi=400)
# %%
# df5.assign(mut_share=lambda x: mut_share(x.day, *popt_wt)).mut_share.plot()
# %%
print('test')
# %%
# %%
xdata = np.linspace(0, 4, 50)


def func(x, a, b, c):
    return a * np.exp(-(b+np.sin(x/360*math.pi)) * x) + c


y = func(xdata, 2.5, 1.3, 0.5)
np.random.seed(1729)
y_noise = 0.2 * np.random.normal(size=xdata.size)
ydata = y + y_noise
popt, pcov = curve_fit(func, xdata, ydata)
popt
# %%
