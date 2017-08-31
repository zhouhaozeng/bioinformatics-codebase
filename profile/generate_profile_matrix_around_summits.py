__author__ = 'Zhouhao Zeng'

import HTSeq
import sys
from optparse import OptionParser
import numpy
import collections


def getSummitProfile(ga, summit_pos_set, window_size, resolution, upstream_extension, downstream_extension):
    upstream_num_points = upstream_extension / resolution
    downstream_num_points = downstream_extension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    profile = collections.defaultdict(lambda: numpy.zeros(total_num_points, numpy.float64))

    for summit_name, summit_pos in summit_pos_set:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            index_pos = summit_pos.pos + (index - upstream_num_points) * resolution
            index_pos_window_iv = HTSeq.GenomicInterval(summit_pos.chrom, index_pos - window_size / 2,
                                                        index_pos + window_size / 2)
            for step_iv, step_count in ga[index_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[summit_name][index] += count_in_window
            index += 1

    return profile


def getSummitProfileWithStrand(ga, summit_pos_set, window_size, resolution, upstream_extension, downstream_extension):
    upstream_num_points = upstream_extension / resolution
    downstream_num_points = downstream_extension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    profile = collections.defaultdict(lambda: numpy.zeros(total_num_points, numpy.float64))

    for summit_name, summit_pos in summit_pos_set:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            if summit_pos.strand == "+":
                index_pos = summit_pos.pos + (index - upstream_num_points) * resolution
                index_pos_window_iv = HTSeq.GenomicInterval(summit_pos.chrom, index_pos - window_size / 2,
                                                            index_pos + window_size / 2)
            elif summit_pos.strand == "-":
                index_pos = summit_pos.pos - (index - upstream_num_points) * resolution
                index_pos_window_iv = HTSeq.GenomicInterval(summit_pos.chrom, index_pos - window_size / 2 + 1,
                                                            index_pos + window_size / 2 + 1)

            for step_iv, step_count in ga[index_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[summit_name][index] += count_in_window
            index += 1

    return profile


def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--tags_bed_file", action="store", type="string", dest="tagsfile", metavar="<file>",
                      help="input ChIP-seq tags bed file")
    parser.add_option("-s", "--summits_file", action="store", type="string", dest="summitsfile", metavar="<file>",
                      help="summits bed file")
    parser.add_option("-t", "--sites_strand_specific", action="store", type="string", dest="sites_strand_specific",
                      metavar="<str>", help="sites strand specific: yes, no", default='no')
    parser.add_option("-n", "--normalization", action="store", type="float", dest="norm", metavar="<float>",
                      help="additional normalization in addition to number of sites, number of reads per million and window_size per 1K")
    parser.add_option("-f", "--fragment_size", action="store", type="int", dest="fragment_size", metavar="<int>",
                      help="fragment size determines the shift (half of fragment_size of ChIP-seq read position, in bps")
    parser.add_option("-r", "--resolution", action="store", type="int", dest="resolution", metavar="<int>",
                      help="resolution of the upstream and downstream profile, eg, 5")
    parser.add_option("-u", "--upstream_extension", action="store", type="int", dest="upstream_extension",
                      metavar="<int>", help="upstream extension")
    parser.add_option("-d", "--downstream_extension", action="store", type="int", dest="downstream_extension",
                      metavar="<int>", help="downstream extension")
    parser.add_option("-w", "--window_size", action="store", type="int", dest="window_size", metavar="<int>",
                      help="window size for averaging. When window size > resolution, there is smoothing")
    parser.add_option("-o", "--output_file", action="store", type="string", dest="outfile", metavar="<file>",
                      help="output profile around sites file")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 20:
        parser.print_help()
        sys.exit(1)

    fragment_size = opt.fragment_size
    window_size = opt.window_size
    resolution = opt.resolution
    upstream_extension = opt.upstream_extension
    downstream_extension = opt.downstream_extension

    print "upstream extension: %i" % upstream_extension
    print "downstream extension: %i" % downstream_extension
    print "upstream and downstream resolution: %i" % resolution
    print "upstream and downstream scanning window size: %i" % window_size

    num_tags = 0
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    tagsfile = HTSeq.BED_Reader(opt.tagsfile)
    for alt in tagsfile:
        if alt.iv.strand == "+":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d + fragment_size / 2)
        elif alt.iv.strand == "-":
            alt_pos = HTSeq.GenomicPosition(alt.iv.chrom, alt.iv.start_d - fragment_size / 2)
        ga[alt_pos] += 1
        num_tags += 1

    summit_pos_set = set()
    summitsfile = HTSeq.BED_Reader(opt.summitsfile)
    for alt in summitsfile:
        summit_pos_set.add((alt.name, alt.iv.start_d_as_pos))

    if opt.sites_strand_specific == "no":
        profile = getSummitProfile(ga, summit_pos_set, window_size, resolution, upstream_extension,
                                   downstream_extension)

    elif opt.sites_strand_specific == "yes":
        profile = getSummitProfileWithStrand(ga, summit_pos_set, window_size, resolution, upstream_extension,
                                             downstream_extension)

    xValues = numpy.arange(-upstream_extension, downstream_extension + 1, resolution)

    f = open(opt.outfile, "w")
    header = 'summit_id' + '\t' + '\t'.join([str(v) for v in xValues]) + '\n'
    f.write(header)

    normalization = num_tags / 1000000.0
    normalization *= window_size / 1000.0
    normalization *= opt.norm

    for summit_name, summit_pos in summit_pos_set:
        profile[summit_name] = profile[summit_name] / normalization

        outline = summit_name + '\t' + '\t'.join([str(v) for v in profile[summit_name]]) + '\n'
        f.write(outline)

    f.close()


if __name__ == "__main__":
    main(sys.argv)
