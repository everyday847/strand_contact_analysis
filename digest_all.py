import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

def make_figs(df, prefix):

    counts = [[0 for aa in aas] for aa in aas]

    # Counts YW and WY separately, but add up to remove this arbitrary order-independence.
    # That means that each WW counts twice. We could alternatively divide everything by two
    for ii,  aa1 in enumerate(aas):
        for jj, aa2 in enumerate(aas):
            if ii > jj: continue
            counts[ii][jj] = len(df[df.arop_aa == f'{aa1}{aa2}']) + len(df[df.arop_aa == f'{aa2}{aa1}'])
            # print(aa1, aa2, counts[ii][jj])

    for ii,  aa1 in enumerate(aas):
        for jj, aa2 in enumerate(aas):
            if ii > jj: counts[ii][jj] = counts[jj][ii]
    
    with open("{}overall_frequencies.csv".format(prefix),'w') as f:
        f.write(","+','.join(aas)+'\n')
        for ii,aa in enumerate(aas):
            f.write(aa+','+','.join([str(counts[ii][jj]) for jj, _aa2 in enumerate(aas)])+'\n')

    ax = sns.heatmap(counts)
    ax.set_xticklabels(aas)
    ax.set_yticklabels(aas)
    plt.savefig('{}overall_frequencies.png'.format(prefix), dpi=300)
    plt.clf()

    # from Parallel and Antiparallel β-Strands Differ in Amino Acid Composition and Availability of Short Constituent Sequences	
    underlying_ap = {
        'A': 6.45,'C':1.63,
        'D':3.26, 'E':5.09, 'F':5.68, 'G':5.19, 'H':2.34, 'I':8.39, 'K':3.12, 'L':9.74,
        'M':1.95, 'N':2.62, 'P':2.14, 'Q':3.16, 'R':5.06, 'S':5.25, 'T':7.02, 'V':12.50, 'W':2.09, 'Y':5.30
    }

    # divide by prevalence
    normalized_counts = {}
    normalized_counts = [[c for c in cc] for cc in counts]
    for idx1 in range(20):
        for idx2 in range(20):
            normalized_counts[idx1][idx2] /= (underlying_ap[aas[idx1]] * underlying_ap[aas[idx2]])

    with open("{}normalized_frequencies.csv".format(prefix),'w') as f:
        f.write(","+','.join(aas)+'\n')
        for ii,aa in enumerate(aas):
            f.write(aa+','+','.join([str(normalized_counts[ii][jj]) for jj, _aa2 in enumerate(aas)])+'\n')

    ax = sns.heatmap(normalized_counts)
    ax.set_xticklabels(aas)
    ax.set_yticklabels(aas)
    plt.savefig('{}normalized_frequencies.png'.format(prefix), dpi=300)
    plt.clf()


if __name__ == '__main__':

    aas = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    df = pd.read_csv('all_pairs_including_nonaromatic.csv')

    # overall -- no redundancy filter
    make_figs(df, '')

    # nonredundant analysis
    # < 2.5A, protein, 30% redundancy
    ids = pd.read_csv('PDBe_search.csv').pdb_id.values

    # filter original df
    df = df[df.pdb.isin(ids)].copy()
    make_figs(df, 'nr_')

    # only hbonded
    dfhb = df[df.arop_is_hbonded == True].copy()
    make_figs(dfhb, 'nr_hb_')

    # only nonhbonded
    dfnhb = df[df.arop_is_hbonded == False].copy()
    make_figs(dfnhb, 'nr_nhb_')


