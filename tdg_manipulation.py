import pandas as pd
import numpy as np

def tdg_virtual4C(tdg, genomic_loci, flank=500000):
    """
    tdg is a pd.DataFrame with columns "x", "y", "z" and index "chr", "pos",
    genomic_loci: a bed like string, e.g. "chr1:100-200"
    flank: the number of base pairs to pull up and down stream of the center
    """
    import re
    chrom, start, end = re.split("[:-]", genomic_loci)
    center = int((int(start) + int(end))/2)
    center_bin = min(tdg.index.get_level_values("pos"), key=lambda x:abs(x-center)) 
    up_stream = center_bin - flank
    down_stream = center_bin + flank
    center_bin_x = tdg.loc[(tdg.index.get_level_values('chr') == chrom) & (tdg.index.get_level_values('pos') == center_bin), "x"].values[0]
    center_bin_y = tdg.loc[(tdg.index.get_level_values('chr') == chrom) & (tdg.index.get_level_values('pos') == center_bin), "y"].values[0]
    center_bin_z = tdg.loc[(tdg.index.get_level_values('chr') == chrom) & (tdg.index.get_level_values('pos') == center_bin), "z"].values[0]

    tdg_subset = tdg.loc[(tdg.index.get_level_values('chr') == chrom) & (tdg.index.get_level_values('pos') >= up_stream) & (tdg.index.get_level_values('pos') <= down_stream)]

    tdg_subset = tdg_subset.assign(distance = lambda x: ((x["x"] - center_bin_x)**2 + (x["y"] - center_bin_y)**2 + (x["z"] - center_bin_z)**2)**0.5)
    return tdg_subset

def distance_between_region(tdg,genomic_loci1,genomic_loci2):
    import re
    chrom1, start1, end1= re.split("[:-]", genomic_loci1)
    center1 = int((int(start1) + int(end1))/2)
    center_bin1 = min(tdg.index.get_level_values("pos"), key=lambda x:abs(x-center1)) 
    center_bin_x1 = tdg.loc[(tdg.index.get_level_values('chr') == chrom1) & (tdg.index.get_level_values('pos') == center_bin1), "x"].values[0]
    center_bin_y1 = tdg.loc[(tdg.index.get_level_values('chr') == chrom1) & (tdg.index.get_level_values('pos') == center_bin1), "y"].values[0]
    center_bin_z1 = tdg.loc[(tdg.index.get_level_values('chr') == chrom1) & (tdg.index.get_level_values('pos') == center_bin1), "z"].values[0]

    chrom2, start2, end2= re.split("[:-]", genomic_loci2)
    center2 = int((int(start2) + int(end2))/2)
    center_bin2 = min(tdg.index.get_level_values("pos"), key=lambda x:abs(x-center2))
    center_bin_x2 = tdg.loc[(tdg.index.get_level_values('chr') == chrom2) & (tdg.index.get_level_values('pos') == center_bin2), "x"].values[0]
    center_bin_y2 = tdg.loc[(tdg.index.get_level_values('chr') == chrom2) & (tdg.index.get_level_values('pos') == center_bin2), "y"].values[0]
    center_bin_z2 = tdg.loc[(tdg.index.get_level_values('chr') == chrom2) & (tdg.index.get_level_values('pos') == center_bin2), "z"].values[0]



    return np.sqrt((center_bin_x1 - center_bin_x2)**2 + (center_bin_y1 - center_bin_y2)**2 + (center_bin_z1 - center_bin_z2)**2