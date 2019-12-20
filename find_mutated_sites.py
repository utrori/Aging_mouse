import collections
import matplotlib.pyplot as plt
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

with open('/mnt/data/old_mouse/rDNA_index/mouse_rDNA_ATP5b2.fasta') as f:
    ref_data = f.readlines()
flag = 0
ref_seq_dict = {}
for line in ref_data:
    if flag == 1:
        ref_seq_dict[refname] = line.strip()
        flag = 0
    elif '>' in line:
        flag = 1
        refname = line[1:].strip()


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


def find_major_mutations(sample, gene, coord1, coord2):
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


def plot_mutations(sample, gene, coord1, coord2):
    """Plot mutation rate.

    Args:
        sample (class): mouse sample class
        gene (str): gene name
        coord1 (int): coordinate
        coord2 (int): coordinate
    """
    bfs = []
    for pos in range(coord1, coord2+1):
        bfs.append(sample.base_freq(gene, pos))

if __name__ == '__main__':
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

    for sample in samplenames:
        print(sample)
        print(len(find_major_mutations(base_summary[sample], '28S', sta_28S, end_28S)))
        plot_mutations(base_summary[sample], '28S', sta_28S, end_28S)
