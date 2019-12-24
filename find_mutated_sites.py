import collections
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob

len_28S = 5491
sta_28S = 177
end_28S = 4906

len_18S = 3570
sta_18S = 263
end_18S = 2132
sta_58S = 3133
end_58S = 3289

len_ATP = 732

gene2len = {'28S': len_28S, '18S+5.8S': len_18S, 'ATP5b': len_ATP}

with open('/mnt/data/old_mouse/rDNA_index/mouse_rDNA_ATP5b2.fasta') as f:
    ref_data = f.readlines()
ref_seq_dict = {}
flag = 0
for line in ref_data:
    if flag:
        ref_seq_dict[refname] = line.strip().upper()
        flag = 0
    elif '>' in line:
        flag = 1
        refname = line[1:].strip()


def determine_consensus_sequence(base_summary):
    """Find consensus sequence of given sample data.

    Args:
        base_summary (dict): a dict from sample name to sample class
    Returns:
        gene2cons: a dictionary that returns consensus sequence from gene name
    """
    genes2cons = {}
    for gene in ref_seq_dict.keys():
        consensus = ''
        for pos in range(gene2len[gene]):
            max_freq_bases = []
            test = []
            for sample in base_summary.values():
                bf = sample.base_freq(gene, pos+1)
                for base, freq in bf.items():
                    if freq > 0.50:
                        max_freq_bases.append(base)
                        test.append((base, freq))
            if len(max_freq_bases) == 6:
                if len(set(max_freq_bases)) == 1:
                    consensus += max_freq_bases[0]
                else:
                    consensus += ref_seq_dict[gene][pos]
            else:
                consensus += ref_seq_dict[gene][pos]
        genes2cons[gene] = consensus
    return genes2cons


def base_repertoire(base_freq_dict):
    """Return which bases are seen as variation.

    Args:
        base_count_dict (dict): return of base_count.
    Returns:
        base_set: a set (literally) of bases seen at the position.
    """
    base_rep = []
    for base, freq in base_freq_dict.items():
        if freq > 0.01:
            base_rep.append(base)
    return set(base_rep)


class mouse_sample():
    """Each sample's base count and frequency are automatically summarized.

    You can access base count and frequency through functions.
    """

    def __init__(self, countdata):
        """Make base count dictionary.

        Args:
            countdata (str): output of bam-readcount
        """
        self.base_count_dict = collections.defaultdict(dict)
        for line in countdata:
            row = line.split()
            gene = row[0]
            pos = int(row[1])
            A_freq = int(line.split('A:')[1].split(':')[0])
            T_freq = int(line.split('T:')[1].split(':')[0])
            G_freq = int(line.split('G:')[1].split(':')[0])
            C_freq = int(line.split('C:')[1].split(':')[0])
            freq_dict = {'A': A_freq, 'T': T_freq,
                         'G': G_freq, 'C': C_freq}
            self.base_count_dict[gene][pos] = freq_dict

    def base_count(self, gene, pos):
        """Return base count.

        This function is added in order keep consistency between frequency and
        count.

        Args:
            gene (str): '28S', '18S+5.8S' or 'ATP5b'
            pos (int): position
        Returns:
            base count in dictionary form
        """
        return self.base_count_dict[gene][pos]

    def base_freq(self, gene, pos):
        """Return base count.

        Args:
            gene (str): '28S', '18S+5.8S' or 'ATP5b'
            pos (int): position
        Returns:
            base frequency in dictionary form
        """
        count_sum = 0
        for count in self.base_count(gene, pos).values():
            count_sum += count
        freq = {}
        for base, count in self.base_count(gene, pos).items():
            freq[base] = count / count_sum
        return freq

    def variation_analysis(self, samples, ref_dict, gene, pos):
        most_freq_base = 'N'
        ref_base = ref_dict[gene][pos-1]
        bf = self.base_freq(gene, pos)
        non_ref_ratio = 1 - bf[ref_base]
        if non_ref_ratio < 0.02:
            return (non_ref_ratio, 'ref')
        base_rep = base_repertoire(bf)
        # base_rep is the repertoire of the sample in question.
        # base_reps is the reportoire of tht other samples
        base_reps = []
        for sample in samples.values():
            base_reps.append(base_repertoire(sample.base_freq(gene, pos)))
        if all(i == base_rep for i in base_reps):
            return (non_ref_ratio, 'common')
        elif all(i == set(ref_base) for i in base_reps) and \
        len(base_rep) > 1:
            return (non_ref_ratio, 'unique')
        else:
            return (non_ref_ratio, 'not_unique')


def minor_mut_rate(sample, gene, coord1, coord2):
    """Return rate of minor mutations in a area of given mouse sample.

    Here, minor mutations refer to mutations at positions,
    where the dominant base is apparent.

    Args:
        sample (class): mouse sample
        gene (str): '28S', '18S+5.8S' or 'ATP5b'
        coord1 (int): start of the coordinate
        coord2 (int): end of the coordinate
    Returns:
        mut_rate: the rate of minor mutations.
    """
    count = 0
    bf_sum = 0
    for pos in range(coord1, coord2 + 1):
        bf = sample.base_freq(gene, pos)
        max_bf = max((i for i in bf.values()))
        if max_bf > 0.9:
            count += 1
            bf_sum += max_bf
    mut_rate = 1 - bf_sum / count
    return mut_rate


def find_major_variations(sample, gene, coord1, coord2):
    """Return a list of major mutations in a area of given mouse sample.

    Here, minor mutations refer to mutations at positions,
    where the dominant base is apparent.

    Args:
        sample (class): mouse sample
        gene (str): '28S', '18S+5.8S' or 'ATP5b'
        coord1 (int): start of the coordinate
        coord2 (int): end of the coordinate
    Returns:
        mut: list of (position, base_frequencey) of mutant positions
    """
    muts = []
    for pos in range(coord1, coord2 + 1):
        bf = sample.base_freq(gene, pos)
        max_bf = max((i for i in bf.values()))
        if max_bf < 0.9:
            muts.append((pos, bf))
    return muts


def analyze_variations(base_summary, ref_dict):
    """Analyze if the mutation is shared by others etc.

    Args:
        base_summary (dict): a dict from sample name to sample class
        ref_dict(dict): dictionary from gene name to the reference sequence
    Returns:
    """
    for gene in ref_seq_dict.keys():
        for pos in range(gene2len[gene]):
            bfs = []
            for sample in base_summary.values():
                bfs.append(sample.base_freq(gene, pos+1))


def plot_mutations(base_summary, sample_name, gene, ref_dict, coord1, coord2):
    """Plot mutation rate.

    Args:
        base_summary (dict): dictionary from sample name to sample class
        sample_name (str): sample name
        gene (str): gene name
        ref_dict(dict): dictionary from gene name to the reference sequence
        coord1 (int): coordinate
        coord2 (int): coordinate
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)

    sample = base_summary[sample_name]
    others = {}
    for k, v in base_summary.items():
        if k != sample_name:
            others[k] = v
    for pos in range(coord1, coord2 + 1):
        value, er_type = sample.variation_analysis(others, ref_dict, gene, pos)
        if er_type == 'ref':
            continue
            ax.plot(pos, value, marker='.', color='black')
        elif er_type == 'common':
            ax.plot(pos, value, marker='.', color=cm.hsv(0.3))
        elif er_type == 'not_unique':
            ax.plot(pos, value, marker='.', color=cm.hsv(0.7))
        elif er_type == 'unique':
            ax.plot(pos, value, marker='.', color=cm.hsv(1))

    ax.set_xlim(coord1, coord2)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Coordinate')
    ax.set_ylabel('Error rate')
    plt.savefig('figs/error_rate_plot_ATP_' + sample_name + '.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    base_summary = pd.read_pickle('base_summary.pkl')
    genes2cons = determine_consensus_sequence(base_summary)
    for sample_name in base_summary.keys():
        plot_mutations(base_summary, sample_name, 'ATP5b', genes2cons, 1, len_ATP)
    quit()
    path = '/mnt/data/old_mouse/*_counted.txt'
    files = glob.glob(path)
    base_summary = {}
    samplenames = []
    for filename in sorted(files):
        samplename = filename.split('/')[-1].split('_')[0]
        samplenames.append(samplename)
        with open(filename) as f:
            countdata = f.readlines()
            base_summary[samplename] = mouse_sample(countdata)

    pd.to_pickle(base_summary, 'base_summary.pkl')
    for sample in samplenames:
        print(sample)
        print(len(find_major_variations(base_summary[sample], '28S', sta_28S, end_28S)))
        plot_mutations(base_summary[sample], '28S', sta_28S, end_28S)
