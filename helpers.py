
def group_tags(taglist):
    """
    For a list of tags like

        ('A', 1), ('A', 2), ...
    
    produce
    
        'A:1-2', ...
    """

    final_tags = []
    chain, start_num, end_num = None, None, None
    seq = ''
    lastchain = None
    for tag in taglist:
        chain, aa, num = tag
        if chain != lastchain:
            final_tags.append((seq, '{}:{}-{}'.format(chain, start_num, end_num)))
            lastchain = chain
            start_num = num
            seq = aa
            end_num = num
            continue
        if num != end_num + 1:
            final_tags.append((seq, '{}:{}-{}'.format(chain, start_num, end_num)))
            lastchain = chain
            start_num = num
            seq = aa
            end_num = num
            continue
        if chain is None:
            chain = tag[0]
        end_num = num
        seq += aa
    
    final_tags.append((seq, '{}:{}-{}'.format(chain, start_num, end_num)))
    return final_tags[1:]

    
def length(strand):
    if '--' not in strand[1] and strand[1].split(':')[1][0] == '-':
        # negative to positive
        nums = strand[1].split(':')[1]
        try:
            return int(nums.split('-')[2]) - (-1 * int(nums.split('-')[1])) + 1
        except:
            print(strand, nums.split('-'))
            quit()
    if '--' in strand[1]:
        try:
            return (-1 * int(strand[1].split(':')[1].split('--')[1])) - int(strand[1].split(':')[1].split('--')[0]) + 1
        except:
            print(strand)
            quit()
    else:
        try:
            return int(strand[1].split(':')[1].split('-')[1]) - int(strand[1].split(':')[1].split('-')[0]) + 1
        except:
            print(strand)
            quit()

            
def find_strands(dssp_fn):
    """
    All strings of E residues are strands.
    """
    strand_lines = []
    with open(dssp_fn) as f:
        flag = False
        for line in f:
            # if line[:len("  #  RESIDUE")+1] == "  #  RESIDUE":
            if "#  RESIDUE" in line:
                flag = True
                continue
            if not flag: continue
            dssr_assignment = line[16]
            chain, aa, seqpos = line[11], line[13], line[5:10].strip()
            if dssr_assignment == 'E':
                try:
                    strand_lines.append((chain, aa, int(seqpos)))
                except ValueError:
                    # print(line)
                    quit()
    
    # They're going to be in order. Simple mapping
    return group_tags(strand_lines)

import gzip
def get_pdb_coords(fn):
    # print(fn)
    
    def chain(line): return line[21]
    def seqpos(line): return int(line[22:26])
    def atomname(line): return line[12:16].strip()
    def coord(line): return (float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()))

    # chain => seqpos => atomname => coords
    coords = {}

    with gzip.open(fn, 'r') as f:
        for line in f:
            line = line.decode('utf8')
            # print(line)
            if line[:6] in ['ATOM  ', 'HETATM']:
                # print(line.strip())
                if chain(line) not in coords:
                    coords[chain(line)] = {}
                if seqpos(line) not in coords[chain(line)]:
                    coords[chain(line)][seqpos(line)] = {}
                coords[chain(line)][seqpos(line)][atomname(line)] = coord(line)
                # print(coords[chain(line)][seqpos(line)][atomname(line)])
    return coords