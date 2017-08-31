__author__ = 'Zhouhao Zeng'

import HTSeq
import sys
from optparse import OptionParser
import numpy


def getSummitProfile(ga, summit_pos_set, window_size, resolution, UpstreamExtension, DownstreamExtension):
    upstream_num_points = UpstreamExtension / resolution
    downstream_num_points = DownstreamExtension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1

    profile = numpy.zeros(total_num_points)
    num_summits = 0

    for summit_pos in summit_pos_set:
        if summit_pos.pos - upstream_num_points * resolution - window_size / 2 < 0:
            continue
        num_summits += 1
        index = 0
        while index < total_num_points:
            count_in_window = 0
            index_pos = summit_pos.pos + (index - upstream_num_points) * resolution
            summit_pos_window_iv = HTSeq.GenomicInterval(summit_pos.chrom, index_pos - window_size / 2,
                                                         index_pos + window_size / 2)

            for step_iv, step_count in ga[summit_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[index] += count_in_window
            index += 1

    return profile, num_summits


def main(argv):
    desc = """This is a template for the analysis of aggretated tag distribution with respect to a set of points, such as the TSSs of known genes, with one profile from each strand."""
    parser = OptionParser(description=desc)
    parser.add_option("-s", "--summits_file", action="store", type="string",
                      dest="summitsfile", metavar="<file>", help="summits bed file")
    parser.add_option("-b", "--bedfile", action="store", type="string",
                      dest="bed_file", help="bed file with interval and score", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="outfile name", metavar="<file>")
    parser.add_option("-u", "--UpstreamExtension", action="store", type="int",
                      dest="upstreamExtension", help="UpstreamExtension", metavar="<int>")
    parser.add_option("-d", "--DownstreamExtension", action="store", type="int",
                      dest="downstreamExtension", help="DownstreamExtension", metavar="<int>")
    parser.add_option("-w", "--WindowSize", action="store", type="int",
                      dest="window_size",
                      help="window size for averaging. When window size > resolution, there is smoothing",
                      metavar="<int>")
    parser.add_option("-r", "--resolution", action="store", type="int",
                      dest="resolution", help="resolution of the upstream and downstream profile, eg, 5",
                      metavar="<int>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 14:
        parser.print_help()
        sys.exit(1)

    window_size = opt.window_size
    UpstreamExtension = opt.upstreamExtension
    DownstreamExtension = opt.downstreamExtension
    resolution = opt.resolution

    print "Upstream extension: %i" % UpstreamExtension
    print "Downstream extension: %i" % DownstreamExtension
    print "Scanning window size: %i" % window_size
    print "Scanning resolution: %i" % resolution

    num_tags = 0
    ga = HTSeq.GenomicArray("auto", stranded=False, typecode="d")
    bedfile = HTSeq.BED_Reader(opt.bed_file)
    for alt in bedfile:
        ga[alt.iv] += alt.score * 1.0 / alt.iv.length

    summit_pos_set = set()
    summitsfile = HTSeq.BED_Reader(opt.summitsfile)
    for alt in summitsfile:
        summit_pos_set.add(alt.iv.start_d_as_pos)

    profile, num_summits = getSummitProfile(ga, summit_pos_set, window_size, resolution, UpstreamExtension,
                                            DownstreamExtension)

    normalization = num_summits
    normalization *= window_size / 1000.0

    print "Number of locations: %i" % num_summits
    print "Normalization = %f" % normalization

    f = open(opt.outfile, "w")
    xValues = numpy.arange(-UpstreamExtension, DownstreamExtension + 1, resolution)
    normalized_profile = profile / normalization
    for index in range(len(xValues)):
        outline = str(xValues[index]) + "\t" + str(normalized_profile[index]) + "\n"
        f.write(outline)
    f.close()


if __name__ == "__main__":
    main(sys.argv)
