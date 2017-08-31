import HTSeq
import sys
from optparse import OptionParser
import numpy
import collections


def getSiteBodyProfile(ga, site_iv_set, genic_partition):
    profile = collections.defaultdict(lambda: numpy.zeros(genic_partition, numpy.float64))
    for site_name, site_iv in site_iv_set:
        partition_size = site_iv.length / genic_partition
        normalization = partition_size * 1.0 / 1000

        index = 0
        for site_pos in site_iv.xrange(partition_size):
            count_in_window = 0
            site_pos_window_iv = HTSeq.GenomicInterval(site_pos.chrom, site_pos.pos, site_pos.pos + partition_size, ".")
            for step_iv, step_count in ga[site_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[site_name][index] += count_in_window / normalization
            index += 1
            if index >= genic_partition:
                break

    return profile


def getUpstreamProfile(ga, site_iv_set, window_size, resolution, upstream_extension, downstream_extension):
    upstream_num_points = upstream_extension / resolution
    downstream_num_points = downstream_extension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    profile = collections.defaultdict(lambda: numpy.zeros(total_num_points, numpy.float64))

    for site_name, site_start_pos in [(site_name, site_iv.start_d_as_pos) for (site_name, site_iv) in site_iv_set]:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            index_pos = site_start_pos.pos + (index - upstream_num_points) * resolution
            index_pos_window_iv = HTSeq.GenomicInterval(site_start_pos.chrom, index_pos - window_size / 2,
                                                        index_pos + window_size / 2)
            for step_iv, step_count in ga[index_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[site_name][index] += count_in_window
            index += 1

    return profile


def getDownstreamProfile(ga, site_iv_set, window_size, resolution, upstream_extension, downstream_extension):
    upstream_num_points = upstream_extension / resolution
    downstream_num_points = downstream_extension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    profile = collections.defaultdict(lambda: numpy.zeros(total_num_points, numpy.float64))

    for site_name, site_end_pos in [(site_name, site_iv.end_d_as_pos) for (site_name, site_iv) in site_iv_set]:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            index_pos = site_end_pos.pos + (index - upstream_num_points) * resolution
            index_pos_window_iv = HTSeq.GenomicInterval(site_end_pos.chrom, index_pos - window_size / 2,
                                                        index_pos + window_size / 2)
            for step_iv, step_count in ga[index_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[site_name][index] += count_in_window
            index += 1

    return profile


def getSiteBodyProfileWithStrand(ga, site_iv_set, genic_partition):
    profile = collections.defaultdict(lambda: numpy.zeros(genic_partition, numpy.float64))
    for site_name, site_iv in site_iv_set:
        partition_size = site_iv.length / genic_partition
        normalization = partition_size * 1.0 / 1000

        index = 0
        for site_pos in site_iv.xrange_d(partition_size):
            count_in_window = 0
            if site_pos.strand == "+":
                site_pos_window_iv = HTSeq.GenomicInterval(site_pos.chrom, site_pos.pos, site_pos.pos + partition_size)
            elif site_pos.strand == "-":
                site_pos_window_iv = HTSeq.GenomicInterval(site_pos.chrom, site_pos.pos - partition_size + 1,
                                                           site_pos.pos + 1)

            for step_iv, step_count in ga[site_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[site_name][index] += count_in_window / normalization
            index += 1
            if index >= genic_partition:
                break

    return profile


def getUpstreamProfileWithStrand(ga, site_iv_set, window_size, resolution, upstream_extension, downstream_extension):
    upstream_num_points = upstream_extension / resolution
    downstream_num_points = downstream_extension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    profile = collections.defaultdict(lambda: numpy.zeros(total_num_points, numpy.float64))

    for site_name, site_start_pos in [(site_name, site_iv.start_d_as_pos) for (site_name, site_iv) in site_iv_set]:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            if site_start_pos.strand == "+":
                index_pos = site_start_pos.pos + (index - upstream_num_points) * resolution
                index_pos_window_iv = HTSeq.GenomicInterval(site_start_pos.chrom, index_pos - window_size / 2,
                                                            index_pos + window_size / 2)
            elif site_start_pos.strand == "-":
                index_pos = site_start_pos.pos - (index - upstream_num_points) * resolution
                index_pos_window_iv = HTSeq.GenomicInterval(site_start_pos.chrom, index_pos - window_size / 2 + 1,
                                                            index_pos + window_size / 2 + 1)

            for step_iv, step_count in ga[index_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[site_name][index] += count_in_window
            index += 1

    return profile


def getDownstreamProfileWithStrand(ga, site_iv_set, window_size, resolution, upstream_extension, downstream_extension):
    upstream_num_points = upstream_extension / resolution
    downstream_num_points = downstream_extension / resolution
    total_num_points = upstream_num_points + downstream_num_points + 1
    profile = collections.defaultdict(lambda: numpy.zeros(total_num_points, numpy.float64))

    for site_name, site_end_pos in [(site_name, site_iv.end_d_as_pos) for (site_name, site_iv) in site_iv_set]:
        index = 0
        while index < total_num_points:
            count_in_window = 0
            if site_end_pos.strand == "+":
                index_pos = site_end_pos.pos + (index - upstream_num_points) * resolution
                index_pos_window_iv = HTSeq.GenomicInterval(site_end_pos.chrom, index_pos - window_size / 2,
                                                            index_pos + window_size / 2)
            elif site_end_pos.strand == "-":
                index_pos = site_end_pos.pos - (index - upstream_num_points) * resolution
                index_pos_window_iv = HTSeq.GenomicInterval(site_end_pos.chrom, index_pos - window_size / 2 + 1,
                                                            index_pos + window_size / 2 + 1)

            for step_iv, step_count in ga[index_pos_window_iv].steps():
                count_in_window += step_count * step_iv.length
            profile[site_name][index] += count_in_window
            index += 1

    return profile


def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--tags_bed_file", action="store", type="string", dest="tagsfile", metavar="<file>",
                      help="input ChIP-seq tags bed file")
    parser.add_option("-s", "--sites_file", action="store", type="string", dest="sitesfile", metavar="<file>",
                      help="sites bed file")
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
    parser.add_option("-p", "--genic_partition", action="store", type="int", dest="genic_partition", metavar="<int>",
                      help="genic partition, eg, 20")
    parser.add_option("-o", "--output_file", action="store", type="string", dest="outfile", metavar="<file>",
                      help="output profile around sites file")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 22:
        parser.print_help()
        sys.exit(1)

    fragment_size = opt.fragment_size
    window_size = opt.window_size
    resolution = opt.resolution
    upstream_extension = opt.upstream_extension
    downstream_extension = opt.downstream_extension
    genic_partition = opt.genic_partition

    print "upstream extension: %i" % upstream_extension
    print "downstream extension: %i" % downstream_extension
    print "upstream and downstream resolution: %i" % resolution
    print "upstream and downstream scanning window size: %i" % window_size
    print "genic partition: %i" % genic_partition

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

    site_iv_set = set()
    sitesfile = HTSeq.BED_Reader(opt.sitesfile)
    for alt in sitesfile:
        site_iv_set.add((alt.name, alt.iv))

    if opt.sites_strand_specific == "no":

        site_body_profile = getSiteBodyProfile(ga, site_iv_set, genic_partition)

        upstream_profile = getUpstreamProfile(ga, site_iv_set, window_size, resolution, upstream_extension, 0)
        upstream_profile *= 1000.0 / window_size

        downstream_profile = getDownstreamProfile(ga, site_iv_set, window_size, resolution, 0, downstream_extension)
        downstream_profile *= 1000.0 / window_size

    elif opt.sites_strand_specific == "yes":

        site_body_profile = getSiteBodyProfileWithStrand(ga, site_iv_set, genic_partition)

        upstream_profile = getUpstreamProfileWithStrand(ga, site_iv_set, window_size, resolution, upstream_extension, 0)
        for site_name, site_iv in site_iv_set:
            upstream_profile[site_name] *= 1000.0 / window_size

        downstream_profile = getDownstreamProfileWithStrand(ga, site_iv_set, window_size, resolution, 0,
                                                            downstream_extension)
        for site_name, site_iv in site_iv_set:
            downstream_profile[site_name] *= 1000.0 / window_size

    upstream_xValues = numpy.arange(0, upstream_extension + 1, resolution)[-1::-1] * (-1)
    site_body_xValues = [0.0] * genic_partition
    for i in xrange(genic_partition):
        site_body_xValues[i] = (i + 0.5) / genic_partition
    downstream_xValues = numpy.arange(0, downstream_extension + 1, resolution)

    f = open(opt.outfile, "w")
    header = 'site_id' + '\t' + '\t'.join([str(v) for v in upstream_xValues]) + '\t' + '\t'.join(
        [str(v) for v in site_body_xValues]) + '\t' + '\t'.join([str(v) for v in downstream_xValues]) + '\n'
    f.write(header)

    normalization = num_tags / 1000000.0
    normalization *= opt.norm

    for site_name, site_iv in site_iv_set:
        site_body_profile[site_name] = site_body_profile[site_name] / normalization
        upstream_profile[site_name] = upstream_profile[site_name] / normalization
        downstream_profile[site_name] = downstream_profile[site_name] / normalization

        outline = site_name + '\t' + '\t'.join([str(v) for v in numpy.hstack(
            [upstream_profile[site_name], site_body_profile[site_name], downstream_profile[site_name]])]) + '\n'
        f.write(outline)

    f.close()


if __name__ == "__main__":
    main(sys.argv)
