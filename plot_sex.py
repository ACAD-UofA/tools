#!/usr/bin/env python

from __future__ import print_function
import sys
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import math
import numpy as np
import operator
from scipy.stats import binom, norm, cauchy, beta
from scipy.cluster.vq import kmeans2


def parse_sextable(filename):
    sextable = []
    with open(filename) as f:
        next(f) # skip header
        for line in f:
            line = line.rstrip()
            fields = line.split("\t")
            if (len(fields) < 5):
                continue
            sample = fields[0]
            Mx = int(fields[1])
            Lx = int(fields[2])
            Ma = int(fields[3])
            La = int(fields[4])

            if Ma < 1000:
                # don't bother
                continue

            Rl = float(Lx)/(Lx+La)
            Rx = float(Mx)/(Mx+Ma)

            # Beta 95%CI with non-informative prior, aka Jefferey's interval.
            # See Brown, Cai, and DasGupta (2001). doi:10.1214/ss/1009213286
            Rx_95CI = beta.interval(0.95, Mx+0.5, Ma+0.5)

            sfields = sample.split("_")
            if sfields[-1].upper() == 'INF':
                age = float('inf')
            else:
                try:
                    age = float(int(sfields[-1]))
                except ValueError:
                    age = None

            sextable.append((sample, Mx, Lx, Ma, La, Rl, Rx, Rx_95CI, age))
    return sextable

def assign_sex(sextable, Rx_init, soft='Cauchy'):
    """assign sex to samples using k-means clustering"""
    Rx = map(operator.itemgetter(6), sextable)
    centroid, hard_classification = kmeans2(Rx, np.array([0.5*Rx_init, Rx_init]), minit='matrix')

    Rx_m = [Rx[i] for i,j in enumerate(hard_classification) if j == 0]
    Rx_f = [Rx[i] for i,j in enumerate(hard_classification) if j == 1]

    m_mu, m_std = centroid[0], np.std(Rx_m)
    f_mu, f_std = centroid[1], np.std(Rx_f)

    if soft is 'Normal':
        # Rx~Normal
        m_dist = norm(m_mu, m_std)
        f_dist = norm(f_mu, f_std)

    elif soft is 'Beta':
        # Rx~Beta - nasty estimation and looks like the normal anyway
        mMx = [Mx[i] for i,j in enumerate(hard_classification) if j == 0]
        mMa = [Ma[i] for i,j in enumerate(hard_classification) if j == 0]
        fMx = [Mx[i] for i,j in enumerate(hard_classification) if j == 1]
        fMa = [Ma[i] for i,j in enumerate(hard_classification) if j == 1]
        mloc, mscale = beta.fit_loc_scale(Rx_m, np.median(mMx), np.median(mMa))
        ma, mb, mloc, mscale = beta.fit(Rx_m, floc=mloc, fscale=mscale)
        #print(ma, mb, mloc, mscale)
        floc, fscale = beta.fit_loc_scale(Rx_f, np.median(fMx), np.median(fMa))
        fa, fb, floc, fscale = beta.fit(Rx_f, floc=floc, fscale=fscale)
        #print(fa, fb, floc, fscale)
        m_dist = beta(ma, mb, loc=mloc, scale=mscale)
        f_dist = beta(fa, fb, loc=floc, scale=fscale)

    elif soft is 'Cauchy':
        # Rx~Cauchy - assumption of independence for num/denom is violated, but otherwise seems sensible
        m_loc, m_scale = cauchy.fit(Rx_m)
        f_loc, f_scale = cauchy.fit(Rx_f)
        m_dist = cauchy(m_loc, m_scale)
        f_dist = cauchy(f_loc, f_scale)
        # use Cauchy central tendency
        m_mu = m_loc
        f_mu = f_loc

    m_bound = m_dist.ppf(0.95)
    f_bound = f_dist.ppf(0.05)
    m_int = m_dist.interval(0.99)
    f_int = f_dist.interval(0.99)
    m_int = (max(0.0,m_int[0]), min(1.0,m_int[1]))
    f_int = (max(0.0,f_int[0]), min(1.0,f_int[1]))

    soft_classification = []
    for x in Rx:
        if x < m_bound:
            soft_classification.append(0)
        elif x > f_bound:
            soft_classification.append(1)
        else:
            soft_classification.append(None)


    if True:
        zn = 100.0
        logit_pm = np.empty_like(Rx)
        cntrd = centroid
        for i, _ in enumerate(Rx):
            Rx_i = Rx[:i] + Rx[i+1:]
            #cntrd, _ = kmeans2(Rx_i, np.array([0.5*Rx_init, Rx_init]), minit='matrix')
            x, a = sextable[i][1], sextable[i][3]
            p_m = 0
            p_f = 0
            for z in np.linspace(m_int[0], m_int[1], zn):
                p_m += binom.pmf(x, x+a, z)*(m_dist.cdf(z+0.5/zn) -m_dist.cdf(z-0.5/zn))
            for z in np.linspace(f_int[0], f_int[1], zn):
                p_f += binom.pmf(x, x+a, z)*(f_dist.cdf(z+0.5/zn) -f_dist.cdf(z-0.5/zn))
            #p_m /= zn
            #p_f /= zn
            if p_m == 0.0:
                logit_pm[i] = -math.log(p_f)
            elif p_f == 0.0:
                logit_pm[i] = math.log(p_m)
            else:
                logit_pm[i] = math.log(p_m) - math.log(p_f)
            #print(p_m, p_f, logit_pm[i])

            #p_m = binom.logsf(x-1, x+a, cntrd[0])
            #p_f = binom.logcdf(x-1, x+a, cntrd[1])
            #logit_pm[i] = p_m - p_f

    return hard_classification, soft_classification, m_dist, f_dist, m_mu, f_mu, m_bound, f_bound, logit_pm

def age_stratification(sextable, sex_assign):
    by_age = filter(lambda z: z[8] is not None, sorted(sextable, key=operator.itemgetter(8)))
    mid_age = int(by_age[1+len(by_age)/2][8])
    young_m = 0
    young_f = 0
    old_m = 0
    old_f = 0
    for i, _a in enumerate(by_age):
        age = _a[8]
        if age is None:
            continue
        if age < mid_age:
            if sex_assign[i] == 0:
                young_m += 1
            elif sex_assign[i] == 1:
                young_f += 1
        else:
            if sex_assign[i] == 0:
                old_m += 1
            elif sex_assign[i] == 1:
                old_f += 1

    # two proportion z-test, normal approx., one sided
    n1 = young_m+young_f
    n2 = old_m+old_f
    p1 = float(young_m)/n1
    p2 = float(old_m)/n2
    phat = (p1*n2 + p2*n1) / (n1 + n2)
    z = (p2 - p1) / math.sqrt(phat*(1.0-phat)*(1.0/n1 + 1.0/n2))
    if z < 0:
        pval = 1.0
    else:
        pval = norm.sf(z)
    #print(n1,n2,p1,p2,phat,z,pval)

    print("{} ({} m, {} f) are younger than {} yBP\n{} ({} m, {} f) are the same age or older than {} yBP".format(young_m+young_f, young_m, young_f, mid_age, old_m+old_f, old_m, old_f, mid_age))
    print("Do older samples have a greater male bias? z={:.2f}, p={:.2f}".format(z, pval))

def print_table(sextable, hard_sex, soft_sex, logit_pm, filename="/dev/stdout"):
    smap = {0:'M', 1:'F', None:'U'}
    with open(filename, "w") as f:
        f.write("sample\tNx\tLx\tNa\tLa\tRx\tRx95CI\tage\tsexhc\tsexsc\tlogitpm\n")
        for _s, hard, soft, lpm in zip(sextable, hard_sex, soft_sex, logit_pm):
            sample, Mx, Lx, Ma, La, _, Rx, Rx95CI, age = _s

            try:
                age = int(age)
            except OverflowError:
                age = "INF"
            except (ValueError, TypeError):
                age = ""

            f.write("{}\t{}\t{}\t{}\t{}\t{:.3f}\t{:.3f}-{:.3f}\t{}\t{}\t{}\t{:.3f}\n".format(
                    sample, Mx, Lx, Ma, La, Rx, Rx95CI[0], Rx95CI[1], age,
                    smap[hard], smap[soft], lpm))

def simulate_sextable(n, sextable, sex_assign):
    from scipy.stats import bernoulli, binom
    import random

    Mx, Lx, Ma, La, Rl, Rx = [np.array(map(operator.itemgetter(i), sextable), dtype=float) for i in range(1,7)]

    Rx_m = np.median([Rx[i] for i,j in enumerate(sex_assign) if j == 0])
    Rx_f = np.median([Rx[i] for i,j in enumerate(sex_assign) if j == 1])

    n_males = len([j for j in sex_assign if j == 0])
    n_females = len([j for j in sex_assign if j == 1])

    # fixed sex ratio
    #alpha = n_males/float(n_males+n_females)
    # ~Beta sex ratio
    alpha = beta.rvs(n_males, n_females)

    rand_table = []
    for i in range(n):
        Nt = int(random.choice(Mx+Ma))
        sex = 1-bernoulli.rvs(alpha)
        if sex == 0:
            # male
            px = beta.rvs(Nt*Rx_m, Nt*(1-Rx_m))
        else:
            # female
            px = beta.rvs(Nt*Rx_f, Nt*(1-Rx_f))
        _Nx = binom.rvs(Nt, px)
        _Na = Nt - _Nx
        _Rx = float(_Nx)/Nt
        _Rx_95CI = beta.ppf([0.025,0.975], _Nx+0.5, _Na+0.5)
        rand_table.append((str(i), _Nx, Lx[0], _Na, La[0], Rl[0], _Rx, _Rx_95CI, None, sex))

    return rand_table

def plot_hist(ax, sextable, mcolour='blue', fcolour='red', ucolour='none'):
    cdict = {0:mcolour, 1:fcolour, None:ucolour}

    Mx, Lx, Ma, La, Rl, Rx = [np.array(map(operator.itemgetter(i), sextable), dtype=float) for i in range(1,7)]
    hard_sex, soft_sex, m_dist, f_dist, m_mu, f_mu, m_bound, f_bound, logit_pm = assign_sex(sextable, sextable[0][5])
    males = [Rx[i] for i,j in enumerate(hard_sex) if j == 0]
    females = [Rx[i] for i,j in enumerate(hard_sex) if j == 1]

    mMx = [Mx[i] for i,j in enumerate(hard_sex) if j == 0]
    mMa = [Ma[i] for i,j in enumerate(hard_sex) if j == 0]
    fMx = [Mx[i] for i,j in enumerate(hard_sex) if j == 1]
    fMa = [Ma[i] for i,j in enumerate(hard_sex) if j == 1]

    hardcolours = [cdict[i] for i in hard_sex]
    softcolours = [cdict[i] for i in soft_sex]
    msoftcolours = [cdict[i] for i in soft_sex if i != 1]
    fsoftcolours = [cdict[i] for i in soft_sex if i != 0]

    _, bins, _ = ax.hist(Rx, bins=len(Rx)/4, fill=False)
    marea = len(males)*(bins[1]-bins[0])
    farea = len(females)*(bins[1]-bins[0])

    x = np.linspace(0, 2*Rl[0], 1000)
    ax.plot(x, marea*m_dist.pdf(x), '-', color=mcolour)
    #ax.scatter(males, marea*m_dist.pdf(males), facecolor=msoftcolours, edgecolor=mcolour)
    ax.plot(x, farea*f_dist.pdf(x), '-', color=fcolour)
    #ax.scatter(females, farea*f_dist.pdf(females), facecolor=fsoftcolours, edgecolor=fcolour)

    ylim = ax.get_ylim()

    ax.vlines(m_mu, 0, ylim[1], color=mcolour, linestyle=':')
    ax.vlines(f_mu, 0, ylim[1], color=fcolour, linestyle=':')
    ax.vlines(m_bound, 0, ylim[1], color=mcolour, linestyle='--')
    ax.vlines(f_bound, 0, ylim[1], color=fcolour, linestyle='--')

    ax.set_ylim((0, ylim[1]))

    m_patch = mpatches.Patch(color=mcolour)
    f_patch = mpatches.Patch(color=fcolour)
    ax.legend([m_patch, f_patch],
            ['male ($n={}({}), \\bar{{x}}={:.3f}, var={:.2e}, Beta.var={:.2e}$)'.format(len(males), len([m for m in males if m<m_bound]), m_mu, np.var(males), beta.var(np.median(mMa), np.median(mMa)+np.median(mMx))),
                'female ($n={}({}), \\bar{{x}}={:.3f}, var={:.2e}, Beta.var={:.2e}$)'.format(len(females), len([f for f in females if f>f_bound]), f_mu, np.var(females), beta.var(np.median(fMa), np.median(fMa)+np.median(fMx)))])
    ax.set_xlabel('Proportion of mapped reads: X/(X+Autosome)')
    ax.set_ylabel('Frequency (counts)')
    ax.set_title('K-means clustering for sex assignment')

def plot_samples(ax, sextable, mcolour='blue', fcolour='red', ucolour='none'):
    cdict = {0:mcolour, 1:fcolour, None:ucolour}
    samples = map(operator.itemgetter(0), sextable)
    Mx, Lx, Ma, La, Rl, Rx = [np.array(map(operator.itemgetter(i), sextable), dtype=float) for i in range(1,7)]
    Rx_95CI = [[_a[7][0] for _a in sextable], [_a[7][1] for _a in sextable]]
    hard_sex, soft_sex, m_dist, f_dist, m_mu, f_mu, m_bound, f_bound, logit_pm = assign_sex(sextable, sextable[0][5])

    ax.vlines(m_mu, -2, len(Rx), color=mcolour, linestyle=':')
    ax.vlines(f_mu, -2, len(Rx), color=fcolour, linestyle=':')
    ax.vlines(m_bound, -2, len(Rx), color=mcolour, linestyle='--')
    ax.vlines(f_bound, -2, len(Rx), color=fcolour, linestyle='--')

    hardcolours = [cdict[i] for i in hard_sex]
    softcolours = [cdict[i] for i in soft_sex]
    msoftcolours = [cdict[i] for i in soft_sex if i != 1]
    fsoftcolours = [cdict[i] for i in soft_sex if i != 0]

    lci = [-lpm if abs(lpm)<10 else -np.sign(lpm)*10 for lpm in logit_pm]

    y_pos = np.arange(len(sextable))
    ax.set_yticks(y_pos)
    ax.set_yticklabels(samples)

    lci = [-lpm if abs(lpm)<10 else -np.sign(lpm)*10 for lpm in logit_pm]
    ax.scatter(Rx, y_pos, facecolor=lci, cmap='bwr', edgecolors='black', lw=0.5, s=60)

    #ax.errorbar(Rx, y_pos, xerr=[Rx-Rx_95CI[0], Rx_95CI[1]-Rx], ecolor=hardcolours, marker="none", fmt="none", capsize=0)

    #ax.annotate("M",
    #        xy=(m_mu, len(sextable)+0.5),
    #        va='top', ha='center',
    #        annotation_clip=False)
    #ax.annotate("F",
    #        xy=(f_mu, len(sextable)+0.5),
    #        va='top', ha='center',
    #        annotation_clip=False)

    ax.set_ylim([-0.5, len(sextable)-0.5])

    ax.set_xlabel('Proportion of mapped reads: X/(X+Autosome)')
    ax.set_ylabel('')
    ax.set_title('Sex assignment')
    #plt.colorbar(sc, ax=ax)


def test_f_std(m, obs_f_std, sextable, hard_sex):

    n_more_extreme = 0
    for i in range(m):
        rtable = simulate_sextable(len(sextable), sextable, hard_sex)
        hs = assign_sex(rtable, Rl[0])[0]
        _Rx = map(operator.itemgetter(6), rtable)
        males = [_Rx[i] for i,j in enumerate(hs) if j == 0]
        females = [_Rx[i] for i,j in enumerate(hs) if j == 1]
        #print(np.std(females)/np.std(males))
        if np.std(females)/np.std(males) >= obs_f_std:
            n_more_extreme += 1

    print(n_more_extreme)
    return n_more_extreme

def test_sex_assign(m, sextable, hard_sex):

    hs_wrong = []
    ss_wrong = []
    ss_wrong2 = []

    for i in range(m):
        rtable = simulate_sextable(len(sextable), sextable, hard_sex)
        hard, soft = assign_sex(rtable, Rl[0])[:2]
        _, soft2 = assign_sex(rtable, Rl[0], soft='Normal')[:2]
        hs_err, ss_err, ss_err2 = 0, 0, 0
        for h, s, s2, _r in zip(hard, soft, soft2, rtable):
            r = _r[9]
            if h != r:
                hs_err += 1
            if s in (0, 1) and s != r:
                ss_err += 1
            if s2 in (0, 1) and s2 != r:
                ss_err2 += 1
        hs_wrong.append(hs_err)
        ss_wrong.append(ss_err)
        ss_wrong2.append(ss_err2)

    from scipy.stats import describe
    print(np.sum(hs_wrong), describe(hs_wrong))
    print('Cauchy', np.sum(ss_wrong), describe(ss_wrong))
    print('Normal', np.sum(ss_wrong2), describe(ss_wrong2))

def plot_regress(ax, x, y, title=None):

    ax.scatter(x, y)

    from scipy.stats import linregress
    slope, intercept, r_val, p_val, std_err = linregress(x, y)
    min_x, max_x = np.min(x), np.max(x)
    min_x = min_x - 0.05*abs(min_x)
    max_x = max_x + 0.05*abs(max_x)

    lx = np.array([min_x, max_x])
    ly = slope*lx + intercept
    ax.plot(lx, ly, '-')

    if title is None:
        title = ""
    else:
        title += ": "
    ax.set_title('{}$r^2$={:.3f}, $p$={:.3f}'.format(title, r_val*r_val, p_val))

    ylim = ax.get_ylim()
    ax.set_ylim([min(0, np.min(y)), ylim[1]])

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} sex.txt out.pdf".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    sextable = parse_sextable(sys.argv[1])
    Mx, Lx, Ma, La, Rl, Rx = [np.array(map(operator.itemgetter(i), sextable), dtype=float) for i in range(1,7)]

    if len(set(Rl)) != 1:
        raise Exception("Bam files not all using the same reference")

    hard_sex, soft_sex, m_dist, f_dist, m_mu, f_mu, m_bound, f_bound, logit_pm = assign_sex(sextable, Rl[0])
    #males = [Rx[i] for i,j in enumerate(hard_sex) if j == 0]
    #females = [Rx[i] for i,j in enumerate(hard_sex) if j == 1]

    #print_table(sextable, hard_sex, soft_sex, logit_pm)
    #exit(1)
    #age_stratification(sextable, soft_sex)

    #from scipy.stats import linregress
    #print(linregress(Rx, Ma))

    #from scipy.stats import mannwhitneyu
    #print(mannwhitneyu(mMx, fMx, alternative='two-sided'))
    #print(mannwhitneyu(mMa, fMa, alternative='two-sided'))


    #test_f_std(10000, np.std(females)/np.std(males), sextable, hard_sex)
    #test_sex_assign(10000, sextable, hard_sex)



    pdf = PdfPages(sys.argv[2])
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(2*fig_w, 6*fig_h))
    gs1 = gridspec.GridSpec(8, 1)
    ax1 = fig1.add_subplot(gs1[0])
    ax2 = fig1.add_subplot(gs1[1],sharex=ax1)
    ax3 = fig1.add_subplot(gs1[2],sharex=ax1)
    ax4 = fig1.add_subplot(gs1[3:],sharex=ax1)

    plot_hist(ax1, sextable)

    plot_regress(ax2, Rx, np.log(Ma), title='Dependency on sequencing depth')
    ax2.set_xlabel('Proportion of mapped reads: X/(X+Autosome)')
    ax2.set_ylabel('log(Autosomal reads)')

    a_Rx, age = zip(*[(st[6], st[8]) for st in sextable if (st[8] is not None and not math.isinf(st[8]))])
    plot_regress(ax3, a_Rx, age, title='Dependency on sample age')
    ax3.set_xlabel('Proportion of mapped reads: X/(X+Autosome)')
    ax3.set_ylabel('Age (years)')

#    for i in range(1,12):
#        ax = fig1.add_subplot(gs1[i])
#        rtable = simulate_sextable(len(sex_table), sextable, hard_sex)
#        plot_hist(ax, rtable)
#        ax.set_xlim(xlim)

    plot_samples(ax4, sextable)

    xlim = [np.min(Rx) -0.05*abs(np.min(Rx)), 1.05*np.max(Rx)]
    ax1.set_xlim(xlim)

    plt.tight_layout()
    pdf.savefig()
    pdf.close()
