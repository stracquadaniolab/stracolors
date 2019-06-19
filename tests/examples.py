import sys
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
import stracolors.settings
import stracolors.calls
import seaborn as sns

# Reset parameters
t=1/10*np.linspace(0,10,100)
x=np.sin(2*np.pi*t)

f,ax=plt.subplots(1)
ax.plot(t,x)
plt.grid(True)
f.savefig('../docs/examples/plot.png', f='png')

stracolors.settings.reset_rc_params()

f,ax=plt.subplots(1)
ax.plot(t,x)
plt.grid(True)
f.savefig('../docs/examples/plot_reset.png', f='png')


# Boxplots

tips = sns.load_dataset("tips")

qual=stracolors.calls.call_palette('qualitative', number_of_colors=4)
f,ax=plt.subplots(1)
sns.boxplot(x="day", y="total_bill",
    data=tips, palette=qual, ax=ax)
f.savefig('../docs/examples/boxplot_qualitative.png', f='png')

paired=stracolors.calls.call_palette_paired(number_of_classes=1)

f,ax=plt.subplots(1)
sns.boxplot(x="day", y="total_bill", hue='smoker',
    data=tips, palette=paired, ax=ax )
f.savefig('../docs/examples/boxplot_paired.png', f='png')


# Triplets of colors
t=t[:np.newaxis]
y1=1+0.8*t
y2=1+t
y3=1+1.2*t

ya=0.2+0.02*np.random.randn(len(t),)
yb=0.5+0.02*np.random.randn(len(t),)
yc=0.8+0.02*np.random.randn(len(t),)

pal=stracolors.calls.call_palette_triple(number_of_classes=2)

f,ax=plt.subplots(1)
plt.plot(t,y1, color=pal[0])
plt.plot(t,y2, color=pal[1])
plt.plot(t, y3, color=pal[2])

plt.plot(t,ya, color=pal[3])
plt.plot(t,yb, color=pal[4])
plt.plot(t, yc, color=pal[5])
f.savefig('../docs/examples/plot_triple.png', f='png')


dots = sns.load_dataset("dots")

# Define a palette to ensure that colors will be
# shared across the facets
palette = dict(zip(dots.coherence.unique(),
                   stracolors.calls.call_palette('sequential',number_of_colors=6)))

# Plot the lines on two facets
g1=sns.relplot(x="time", y="firing_rate",
            hue="coherence", size="choice", col="align",
            size_order=["T1", "T2"], palette=palette,
            height=5, aspect=.75, facet_kws=dict(sharex=False),
            kind="line", legend="full", data=dots)
g1.fig.savefig('../docs/examples/sequential_red.png', f='png')


# Define a palette to ensure that colors will be
# shared across the facets
palette = dict(zip(dots.coherence.unique(),
                   stracolors.calls.call_palette('sequential',number_of_colors=6, shade='blue', reverse=True)))
# Plot the lines on two facets
g1=sns.relplot(x="time", y="firing_rate",
            hue="coherence", size="choice", col="align",
            size_order=["T1", "T2"], palette=palette,
            height=5, aspect=.75, facet_kws=dict(sharex=False),
            kind="line", legend="full", data=dots)
g1.fig.savefig('../docs/examples/sequential_blue.png', f='png')