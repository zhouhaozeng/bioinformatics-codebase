import sys, os, re
from optparse import OptionParser
import subprocess

def main(argv):
    parser = OptionParser()
    parser.add_option("-p", "--promoter_iv_file", action="store", type="string", dest="promoterfile", metavar="<file>", help="promoter iv file")
    parser.add_option("-l", "--histone_mods_list", action="store", type="string", dest="hmlist", metavar="<file>", help="histone mods list")
    parser.add_option("-e", "--enhancer_file_prefix", action="store", type="string", dest="prefix", metavar="<file>", help="Output enhancer files prefix")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 6:
        parser.print_help()
        sys.exit(1)    
    
    hm_file_dict = {}
    with open(opt.hmlist, 'r') as f:
        for line in f:
            line = line.strip()
            sline = line.split()
            hm_file_dict[sline[0]] = sline[1]
    
            
    str_cmd = 'bedtools subtract -a %s -b %s > %s' % (hm_file_dict['H3K4Me1'], hm_file_dict['H3K4Me3'], opt.prefix + '_Me1+Me3-.bed')
    p = subprocess.Popen(str_cmd, shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    
    str_cmd = 'bedtools intersect -v -a %s -b %s > %s' % (opt.prefix + '_Me1+Me3-.bed', opt.promoterfile, opt.prefix + '_enhancer.bed')
    p = subprocess.Popen(str_cmd, shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    
    str_cmd = 'bedtools intersect -v -a %s -b %s | bedtools intersect -v -a stdin -b %s > %s' % (opt.prefix + '_enhancer.bed', hm_file_dict['H3K27Ac'], hm_file_dict['H3K27Me3'], opt.prefix + '_primed_enhancer.bed')
    p = subprocess.Popen(str_cmd, shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    
    str_cmd = 'bedtools intersect -u -a %s -b %s | bedtools intersect -v -a stdin -b %s > %s' % (opt.prefix + '_enhancer.bed', hm_file_dict['H3K27Ac'], hm_file_dict['H3K27Me3'], opt.prefix + '_active_enhancer.bed')
    p = subprocess.Popen(str_cmd, shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    
    str_cmd = 'bedtools intersect -u -a %s -b %s | bedtools intersect -v -a stdin -b %s > %s' % (opt.prefix + '_enhancer.bed', hm_file_dict['H3K27Me3'], hm_file_dict['H3K27Ac'], opt.prefix + '_poised_enhancer.bed')
    p = subprocess.Popen(str_cmd, shell=True, stdout = subprocess.PIPE, stderr= subprocess.PIPE)
    p.wait()
    
if __name__ == "__main__":
    main(sys.argv)
            
            
            
    