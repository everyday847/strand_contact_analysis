#!/usr/bin/env python
from tqdm import tqdm
from helpers import group_tags, length, find_strands, get_pdb_coords
import pandas as pd 
import sys

"""
Finds all antiparallel beta strands in every PDB. 
Finds the ones with *any* amino acid pairs.
"""


def at_least_one_aromatic(strands):
    return [strand for strand in strands if 'F' in strand[0] or 'W' in strand[0] or 'Y' in strand[0]]


def at_least_three_aa(strands):
    return [strand for strand in strands if length(strand) >= 3]




import glob

def dsq(c1, c2):
    return (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2

def termini(strand):
    nums = strand[1].split(':')[1]

    if '--' not in nums and nums[0] == '-':
        return (-1 * int(nums.split('-')[1]), int(nums.split('-')[2]))
    if '--' in nums:
        return (int(nums.split('--')[0]), -1 * int(nums.split('--')[1]))
    else:
        return (int(nums.split('-')[0]), int(nums.split('-')[1]))

def no_pair_CA_close(strand1, strand2, coords):
    """
    Quick filter on strands to ensure that we only consider strand pairs
    with at least one pair of CA atoms reasonably nearby. There are way
    more strand pairs than strand pairs near each other in any given
    structure.
    """

    strand1c = strand1[1][0]
    strand1s, strand1e = termini(strand1)
    strand2c = strand2[1][0]
    strand2s, strand2e = termini(strand2) 
    for s1i in range(strand1s, strand1e+1):
        for s2i in range(strand2s, strand2e+1):
            # # print(strand1c, s1i, strand2c, s2i)
            # print(coords[strand1c][s1i]['CA'],
            #     coords[strand2c][s2i]['CA'])
            try:
                dist_sq = dsq(
                    coords[strand1c][s1i]['CA'],
                    coords[strand2c][s2i]['CA'],
                )
            except KeyError:
                # print(strand1c, s1i, strand2c, s2i, end='')
                raise KeyError
            # print(dist_sq)
            if dist_sq <= 100: return False
            
            # handle NMR structures where there could be overlapping chains due
            # to accidental model concatenation during processing
            if dist_sq < 1: return True

    return True

def h_bonds(strand1, strand2, coords):
    """
    Finds all hydrogen bonds between two strands. Uses a pretty
    generous cutoff to find poorly oriented hbonds.
    """

    hbs = []
    strand1c = strand1[1][0]
    strand1s, strand1e = termini(strand1)
    strand2c = strand2[1][0]
    strand2s, strand2e = termini(strand2)
    for c1, s1i in zip(strand1[0], range(int(strand1s), int(strand1e)+1)):
        for c2, s2i in zip(strand2[0], range(int(strand2s), int(strand2e)+1)):
            # print(coords[strand1c][s1i])
            dist_sq = dsq(
                coords[strand1c][s1i]['N'],
                coords[strand2c][s2i]['O'],
            )
            if dist_sq < 12.25:
                hbs.append('{}:{} and {}:{}'.format(strand1c, s1i, strand2c, s2i))
                # print('possible h-bond between', strand1c, s1i, 'N', strand2c, s2i, 'O')

            dist_sq = dsq(
                coords[strand1c][s1i]['O'],
                coords[strand2c][s2i]['N'],
            )
            if dist_sq < 12.25:
                hbs.append('{}:{} and {}:{}'.format(strand1c, s1i, strand2c, s2i))
                #print('possible h-bond between', strand1c, s1i, 'O', strand2c, s2i, 'N')

    return hbs

                
def aro_pairs(strand1, strand2, coords):
    strand1c = strand1[1][0]
    strand1s, strand1e = termini(strand1)
    strand2c = strand2[1][0]
    strand2s, strand2e = termini(strand2)
    arop = []
    arop_aa = []
    for c1, s1i in zip(strand1[0], range(int(strand1s), int(strand1e)+1)):
        for c2, s2i in zip(strand2[0], range(int(strand2s), int(strand2e)+1)):
            # print(coords[strand1c][s1i])
            # if c1 not in 'FYW' or c2 not in 'FYW': continue
            # if 'CB' not in coords[strand1c][s1i] or 'CB' not in coords[strand2c][s2i]: continue
            # dist_sq = dsq(
            #     coords[strand1c][s1i]['CB'],
            #     coords[strand2c][s2i]['CB'],
            # )
            # if dist_sq < 25:
            #     print('possible sc contact', strand1c, s1i, 'N', strand2c, s2i, 'O')
            done_with_residue_pair = False
            for k1, v1 in coords[strand1c][s1i].items():
                if k1.strip() in ['H', 'N', 'O', 'C']: continue
                for k2, v2 in coords[strand2c][s2i].items():
                    if k2.strip() in ['H', 'N', 'O', 'C']: continue

                    dist_sq = dsq(v1, v2)
                    if dist_sq < 16:
                        has_ap = True
                        # arop_aa.append('{}{}'.format(c1, c2))
                        arop_aa.append('{}{}'.format(min([c1, c2]), max([c1, c2])))
                        arop.append('{}:{} and {}:{}'.format(strand1c, s1i, strand2c, s2i))
                        done_with_residue_pair = True
                        break
                if done_with_residue_pair: break
                        # return True, 

    return arop, arop_aa

def orient(hbs):
    """
    In parallel beta strands, the sequence distance between hydrogen bonded
    residues remains comparable. Maybe there's a bulge here and there, but
    there are no big changes. In contrast, antiparallel beta strands have
    a different sequence distance for each hbond -- one end is as near as
    can be and the other is as far as can be. This easily classifies
    the two orientations.
    """
    def seqdiff(hb):
        rs = list(map(int, map(lambda r: r.split(':')[1], hb.split(' and '))))
        return rs[1] - rs[0]

    sequence_diffs = []
    for hb in hbs:
        sequence_diffs.append(seqdiff(hb))
    # print(sequence_diffs)
    # return 'p' if len(set(sequence_diffs)) == 1 else 'a'
    return 'p' if max(sequence_diffs) - min(sequence_diffs) < 3 else 'a'


import multiprocessing

def data_of(pdb_dir, fns):
    """
    Grab all the salient data from a batch of filenames
    """

    # print(jobid)
    data = []
    for fn in tqdm(fns):
        # print(fn[-9:-5])
        # print(fn)

        # print(at_least_one_aromatic(find_strands('/Users/andrewwatkins/data/dssp/{}.dssp'.format(fn[-11:-7]))))
        valid_strands = at_least_three_aa(find_strands(fn))
        # print(len(valid_strands))
        # look in PDB for coordinates and stuff
        try:
            coords = get_pdb_coords(f'{pdb_dir}/{fn[-8:-6]}/pdb{fn[-9:-5]}.ent.gz')
        except FileNotFoundError: 
            # print('/Users/andrewwatkins/data/pdb/{}/pdb{}.ent.gz'.format(fn[-8:-6], fn[-9:-5]))
            continue
        
        for ii, strand1 in enumerate(valid_strands):
            for jj, strand2 in enumerate(valid_strands[ii+1:]):
                # if strand1[0] == strand2[0] and strand1[1] == strand1[1]: continue
                try:
                    if no_pair_CA_close(strand1, strand2, coords): continue
                except KeyError:
                    # print('({})'.format(fn[-9:-5]))
                    continue
                try:
                    hbs = h_bonds(strand1, strand2, coords)
                except KeyError:
                    # print('({})'.format(fn[-9:-5]))
                    continue
                if len(hbs) == 0: continue
                orientation = orient(hbs)
                if orientation == 'p': continue
                arop, arop_aa = aro_pairs(strand1, strand2, coords)
                # print(has_ap, ap)

                for arop_, arop_aa_ in zip(arop, arop_aa):
                    if arop_ in hbs: arop_is_hbonded = True
                    elif ' and '.join(arop_.split(' and ')[::-1]) in hbs: 
                        # print(arop)
                        arop_is_hbonded = True
                    else:
                        arop_is_hbonded = False
                    
                    data.append({
                        'pdb': fn[-9:-5],
                        'res1': strand1[1],
                        'res2': strand2[1],
                        'seq1': strand1[0],
                        'seq2': strand2[0],
                        'arop': arop_,
                        'arop_aa': arop_aa_,
                        'arop_is_hbonded': arop_is_hbonded
                    })
                # except IndexError:
                #     print(arop)
                #     quit()
    return data

if __name__ == '__main__':
    data = []

    def chunks(l, n):
        return [l[i:i+n] for i in range(0, len(l), n)]
    
    if len(sys.argv) != 3:
        print('Usage: python main_all_contacts.py {dssp_dir} {pdb_dir}')
        print('Pass a directory where you have dssp files and pdb files, respectively.')
        print('Each directory should have subdirectories named after the middle two letters of the PDB code.')

    dssp_dir = sys.argv[1]
    pdb_dir = sys.argv[2]
    print(dssp_dir)

    l = list(glob.glob(f'{dssp_dir}/*.dssp'))

    nproc = 16
    total = len(l)
    chunk_size = total // nproc
    tslice = chunks(l, chunk_size)
    jobs = []
    with multiprocessing.Pool(processes=nproc) as pool:
        data = pool.map(data_of, tslice)
            # for i, s in enumerate(tslice):
    #     j = multiprocessing.Process(target=data_of, args=(i, s))
    #     jobs.append(j)
    # for j in jobs:
    #     j.start()

    # dispatch_jobs(l, 16)

    df = pd.concat([pd.DataFrame.from_dict(d) for d in data])
    df.to_csv('all_pairs_including_nonaromatic.csv')       
