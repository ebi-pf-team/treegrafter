#!/usr/bin/env python

import os
import sys
import re
import json
import argparse
import logging
from Bio import Phylo
from Bio.Phylo import NewickIO

# relative imports
from tglib.re_matcher import re_matcher

# options = []



def _querymsf(matchdata, pthrAlignLength):

    # matchdata contains: hmmstart, hmmend, hmmalign and matchalign, as arrays (multiple modules possible)


    # N-terminaly padd the sequence
    # position 1 until start is filled with '-'
    querymsf = ((int(matchdata['hmmstart'][0]) - 1) * '-')

    # loop the elements/domains
    for i in range(0, len(matchdata['matchalign'])):

        # if this is not the first element, fill in the gap between the hits
        if i > 0:
            start = int(matchdata['hmmstart'][i])
            end = int(matchdata['hmmend'][i-1])
            # This bridges the query_idgap between the hits
            querymsf += (start - end -1) * '-'

        # extract the query string
        matchalign = matchdata['matchalign'][i]
        hmmalign = matchdata['hmmalign'][i]

        # loop the sequence
        for j in range(0, len(hmmalign)):
            # hmm insert state
            if hmmalign[j:j+1] == ".":
                continue

            querymsf += matchalign[j:j+1]

    # C-terminaly padd the sequence
    # get the end of the last element/domain
    last_end = int(matchdata['hmmend'][-1])
    # and padd out to fill the msf lenght
    querymsf += (int(pthrAlignLength) - last_end) * '-'

    # error check (is this required?)
    if (len(querymsf) != pthrAlignLength):
        # then something is wrong

        print(matchdata)
        print(querymsf)
        print(pthrAlignLength)
        print(len(querymsf))

        return 0

    return querymsf.upper()


def stringify(query_id):
    # stringify query_id
    
    query_id = re.sub('[^\w]', '_', query_id)

    return query_id


def _generateFasta(pathr, query_id, querymsf):

    # use static dir paths for testing. these are provided
    fasta_in_dir = options['msf_tree_folder']
    fasta_out_dir = options['tmp_folder']
    fasta_out_file = fasta_out_dir + query_id + '.' + pathr + '.fasta'

    with open(fasta_out_file, 'w') as outfile:
        with open(fasta_in_dir + pathr + '.AN.fasta', 'r') as infile:
            # print first fasta header
            outfile.write('>query_' + query_id + '\n')
            # and the body lines
            for i in range(0, len(querymsf), 80):
                outfile.write(querymsf[i:i+80] + '\n')
            # copy over the other fasta file
            outfile.write(infile.read())

    return(fasta_out_file)


def _run_raxml(pathr, query_id, fasta_file, annotations):

    # print(fasta_file)

    pantherdir = options['msf_tree_folder']

    bifurnewick_in = pantherdir + pathr + '.bifurcate.newick'

    raxml_dir = options['tmp_folder'] + \
        pathr + '_' + query_id + '_raxml' + str(os.getpid())

    os.mkdir(raxml_dir)

    exit_status = os.system('raxmlHPC-PTHREADS-SSE3 -f y -p 12345 -t ' + bifurnewick_in +
                            ' -G 0.05 -m PROTGAMMAWAG -T 4 -s ' + fasta_file + ' -n ' + pathr + ' -w ' + raxml_dir + ' >/dev/null')

    # print(exit_status)

    mapANs = _mapto(raxml_dir, pathr, query_id)
    # print(mapANs)

    commonAN = _commonancestor(pathr, mapANs)
    # print(commonAN)

    if commonAN is None:
        commonAN = 'root'

    annot = annotations[pathr + ':' + str(commonAN)]
    # print(annot)
    result = query_id + "\t" + pathr + "\t" + annot + "\n"
    # print(result)
    return(result)




def _mapto(raxml_dir, pathr, query_id):

    # print(raxml_dir, pathr, query_id)

    classification_file = raxml_dir + '/RAxML_portableTree.' + pathr + '.jplace'
    # print(classification_file)

    with open(classification_file) as classification:
        classification_json = json.load(classification)

    # print(classification_json)


    tree_string = classification_json['tree']
    # print(tree_string)

    # tree_string = "(((((AN11:0.84399999999999997247{0},((AN7:1.00600000000000000533{1},AN8:0.87399999999999999911{2}):0.14999999999999999445{3},(AN9:0.77700000000000002398{4},AN10:2.00000000000000000000{5}):0.14999999999999999445{6}):0.14999999999999999445{7}):0.00500000000000000010{8},(AN13:0.63000000000000000444{9},AN14:0.50200000000000000178{10}):0.04700000000000000011{11}):0.65000000000000002220{12},(AN3:1.52200000000000001954{13},AN4:1.17700000000000004619{14}):0.14999999999999999445{15}):0.24399999999999999467{16},(((AN17:1.04200000000000003730{17},((AN20:0.53600000000000003197{18},AN21:0.55800000000000005151{19}):0.56799999999999994937{20},AN22:0.79700000000000004174{21}):0.46000000000000001998{22}):0.55900000000000005240{23},(AN24:0.70499999999999996003{24},AN25:0.64400000000000001688{25}):0.78500000000000003109{26}):0.14999999999999999445{27},((AN27:0.71899999999999997247{28},AN28:0.58299999999999996270{29}):0.63300000000000000711{30},((AN31:0.71899999999999997247{31},(AN37:0.16400000000000000688{32},((AN34:0.83899999999999996803{33},AN35:0.17299999999999998712{34}):0.05999999999999999778{35},AN36:0.36399999999999999023{36}):0.14999999999999999445{37}):0.54800000000000004263{38}):0.43800000000000000044{39},(((AN45:0.14099999999999998646{40},(AN43:0.17000000000000001221{41},AN44:0.47399999999999997691{42}):0.14999999999999999445{43}):0.58899999999999996803{44},(AN49:0.05899999999999999689{45},(AN47:0.07199999999999999456{46},AN48:0.25500000000000000444{47}):0.14999999999999999445{48}):0.43499999999999999778{49}):0.20799999999999999045{50},(AN39:0.51200000000000001066{51},AN40:0.44500000000000000666{52}):0.14999999999999999445{53}):0.64500000000000001776{54}):0.48099999999999998312{55}):0.14999999999999999445{56}):0.25100000000000000089{57}):0.39800000000000002043{58},((((AN75:0.75400000000000000355{59},(((AN79:0.17799999999999999156{60},((AN82:0.22300000000000000377{61},AN83:0.20999999999999999223{62}):0.52800000000000002487{63},(AN85:0.33100000000000001643{64},AN86:0.30799999999999999600{65}):0.02400000000000000050{66}):0.18399999999999999689{67}):0.12199999999999999734{68},(AN88:0.14799999999999999267{69},(((AN94:0.19600000000000000755{70},AN95:0.06099999999999999867{71}):0.08200000000000000344{72},(AN97:2.00000000000000000000{73},AN98:0.42499999999999998890{74}):0.05500000000000000028{75}):0.40600000000000002753{76},(AN90:0.16300000000000000600{77},AN91:0.17299999999999998712{78}):0.14999999999999999445{79}):0.08300000000000000433{80}):0.13500000000000000888{81}):0.22700000000000000733{82},(AN100:2.00000000000000000000{83},AN101:0.48899999999999999023{84}):0.09600000000000000200{85}):0.38000000000000000444{86}):0.47999999999999998224{87},(((AN55:0.67200000000000004174{88},AN56:0.60799999999999998490{89}):0.58599999999999996536{90},(((AN59:1.10600000000000009415{91},AN60:1.22199999999999997513{92}):0.14999999999999999445{93},(AN61:1.24300000000000010481{94},AN62:1.10600000000000009415{95}):0.14999999999999999445{96}):0.17599999999999998979{97},(AN64:0.17199999999999998623{98},AN65:0.18699999999999999956{99}):0.45800000000000001821{100}):0.66900000000000003908{101}):0.46800000000000002709{102},((AN70:1.05299999999999993605{103},(AN68:1.24700000000000010836{104},AN69:1.20199999999999995737{105}):0.14999999999999999445{106}):0.12099999999999999645{107},(AN72:0.90100000000000002309{108},AN73:0.92200000000000004174{109}):0.13600000000000000977{110}):0.49699999999999999734{111}):0.14999999999999999445{112}):0.29999999999999998890{113},(((AN105:0.76900000000000001688{114},AN106:0.73599999999999998757{115}):0.16200000000000000511{116},(AN110:0.67500000000000004441{117},(AN108:0.55500000000000004885{118},AN109:0.76800000000000001599{119}):0.14999999999999999445{120}):0.32300000000000000933{121}):0.58499999999999996447{122},AN111:1.19999999999999995559{123}):0.39400000000000001688{124}):0.13400000000000000799{125},(((((AN116:1.02800000000000002487{126},AN117:0.98099999999999998312{127}):0.32600000000000001199{128},AN118:0.65400000000000002576{129}):0.21900000000000000022{130},((AN121:0.14499999999999999001{131},AN122:0.14499999999999999001{132}):0.38600000000000000977{133},((AN125:0.07699999999999999900{134},AN126:0.04299999999999999656{135}):0.03300000000000000155{136},AN127:0.02000000000000000042{137}):0.45800000000000001821{138}):0.44300000000000000488{139}):1.09400000000000008349{140},(AN129:1.00000000000000000000{141},(AN131:0.22500000000000000555{142},AN132:0.20399999999999998690{143}):1.01499999999999990230{144}):1.13999999999999990230{145}):0.14999999999999999445{146},((((AN135:0.74099999999999999201{147},((AN138:0.96899999999999997247{148},AN139:0.76400000000000001243{149}):0.25100000000000000089{150},((AN142:0.60199999999999997957{151},AN143:0.57099999999999995204{152}):0.13400000000000000799{153},(AN145:0.58699999999999996625{154},AN146:0.49199999999999999289{155}):0.13800000000000001155{156}):0.31900000000000000577{157}):0.21900000000000000022{158}):0.56000000000000005329{159},(AN148:0.57299999999999995381{160},AN149:0.58099999999999996092{161}):0.33400000000000001910{162}):0.14999999999999999445{163},(AN150:1.11499999999999999112{164},AN151:0.77800000000000002487{165}):0.14999999999999999445{166}):1.09099999999999996980{167},((AN206:0.95499999999999996003{168},((AN213:0.04599999999999999922{169},(AN211:0.03799999999999999906{170},AN212:0.39000000000000001332{171}):0.14999999999999999445{172}):0.54300000000000003819{173},(AN208:0.85999999999999998668{174},AN209:0.51200000000000001066{175}):0.14999999999999999445{176}):0.35399999999999998135{177}):0.46200000000000002176{178},(((((AN157:0.82799999999999995826{179},(AN155:1.17799999999999993605{180},AN156:1.09400000000000008349{181}):0.14999999999999999445{182}):0.08200000000000000344{183},(((AN179:0.88800000000000001155{184},(AN177:1.03600000000000003197{185},AN178:0.98799999999999998934{186}):0.14999999999999999445{187}):0.06500000000000000222{188},(AN181:0.75000000000000000000{189},(AN183:0.62600000000000000089{190},(AN187:1.00499999999999989342{191},(AN185:1.33200000000000007283{192},AN186:1.50699999999999989519{193}):0.14999999999999999445{194}):0.17199999999999998623{195}):0.31800000000000000488{196}):0.11200000000000000233{197}):0.35299999999999998046{198},((AN160:0.80500000000000004885{199},(AN162:0.50900000000000000799{200},(AN164:0.47899999999999998135{201},AN165:0.40799999999999997380{202}):0.25900000000000000799{203}):0.44500000000000000666{204}):0.63200000000000000622{205},(AN167:1.00200000000000000178{206},((AN170:1.21599999999999996980{207},AN171:0.52400000000000002132{208}):0.05500000000000000028{209},(AN173:0.83399999999999996358{210},AN174:0.50600000000000000533{211}):0.05399999999999999939{212}):0.56399999999999994582{213}):0.81200000000000005507{214}):0.14999999999999999445{215}):0.00800000000000000017{216}):0.64600000000000001865{217},(((AN191:0.38400000000000000799{218},AN192:0.34799999999999997602{219}):0.30699999999999999512{220},AN193:0.53500000000000003109{221}):0.31300000000000000044{222},AN194:0.58599999999999996536{223}):0.70399999999999995914{224}):0.14999999999999999445{225},(AN195:1.05499999999999993783{226},(AN197:0.88500000000000000888{227},AN198:0.65200000000000002398{228}):0.69799999999999995381{229}):0.14999999999999999445{230}):0.14999999999999999445{231},((AN199:1.52000000000000001776{232},AN200:1.15999999999999992006{233}):0.14999999999999999445{234},((AN202:0.45700000000000001732{235},AN203:0.48399999999999998579{236}):0.64800000000000002043{237},AN204:1.56699999999999994849{238}):0.14999999999999999445{239}):0.14999999999999999445{240}):0.14999999999999999445{241}):0.85399999999999998135{242}):0.14999999999999999445{243}):0.05500000000000000028{244}):0.39800000000000002043{245});"

    matches = re.findall('AN(\d+):\d+\.\d+\{(\d+)\}', tree_string)
    # print(matches)

    AN_label = {}
    for [an, r] in matches:
        AN_label['AN' + an] = 'R' + r
        AN_label['R' + r] = 'AN' + an

    # print(AN_label)

    newick_string = re.sub('(AN\d+)?\:\d+\.\d+{(\d+)}', 'R\g<2>', tree_string)
    # print(newick_string)


    mytree = Phylo.read(NewickIO.StringIO(newick_string), 'newick')
    # print(mytree)
    # Phylo.draw_ascii(mytree)


    locations_ref = classification_json['placements'][0]['p']
    # locations_ref = [[130, 13902], [238, 13902]]
    # print(locations_ref)

    child_ids = []

    ter = []

    for maploc in locations_ref:
        # print("maploc")

        rloc = 'R' + str(maploc[0])

        # print(rloc)

        node = mytree.find_clades(rloc).__next__()
        # print(node)

        ter.extend(node.get_terminals())

    # print("maploc OUT")

    comonancestor = mytree.common_ancestor(ter)

    # print(comonancestor)

    for leaf in comonancestor.get_terminals():
        child_ids.append(AN_label[leaf.name])

    # print(child_ids)
    return child_ids
 

def _commonancestor(pathr, mapANs):

    pantherdir = options['msf_tree_folder']

    newick_in = pantherdir + pathr + '.newick'

    newtree = Phylo.read(newick_in, "newick")
    # Phylo.draw_ascii(newtree)

    commonancestor = newtree.common_ancestor(mapANs)
    # print(commonancestor)

    return commonancestor


def runhmmr():
    options['hmmr_out'] = options['fasta_input'] + '.' + options['hmmr_mode'] + '.out'

    panther_hmm = options['data_folder'] + \
        '/famhmm/binHmm'

    hmmr_cmd = options['hmmr_mode'] + \
        ' --notextw --cpu ' + \
        str(options['hmmr_cpus']) + ' -o ' + options['hmmr_out'] + \
        ' ' + panther_hmm + ' ' + options['fasta_input'] + ' > /dev/null'

    if options['hmmr_bin'] is not None:
        hmmr_cmd = options['hmmr_bin'] + '/' + hmmr_cmd

    exit_status = os.system(hmmr_cmd)

    if exit_status != 1:
        sys.exit('Error running hmmer')

    return exit_status



def parsehmmr(hmmer_out):
    if options['hmmr_mode'] == 'hmmscan':
        return parsehmmscan(hmmer_out)
    elif options['hmmr_mode'] == 'hmmsearch':
        return parsehmmsearch(hmmer_out)


def parsehmmscan(hmmer_out):

    matches = {}
    

    with open(hmmer_out) as fp:
        align_found_n = 0
        matchpthr = None
        query_id = None

        line = fp.readline()
        while line:
            m = re_matcher(line)
            # print("Line {}: {}".format(cnt, line.strip()))

            if line.startswith('#') or not line.strip():
                # print("INSIDE IF 1: {}".format(line.strip()))
                line = fp.readline()
                continue
            elif m.match('\AQuery:\s+(\S+)'):
                # print("INSIDE IF 2: {}".format(line.strip()))
                query_id = m.group(1)
                # print(query_id)
            elif m.match('\A>> (PTHR[0-9]+)'):
                align_found_n += 1
                # print("INSIDE IF 3: {}".format(line.strip()))
                matchpthr = m.group(1)
                # print(matchpthr)
                # print(align_found_n)
            elif m.match('\s+\d+\s+!') and align_found_n == 1:
                # print("INSIDE IF 4: {}".format(line.strip()))

                mark = line.split()
                # print(mark)

                if matchpthr not in matches:
                    matches[matchpthr] = {}
                if query_id not in matches[matchpthr]:
                    matches[matchpthr][query_id] = {}
                if 'hmmstart' not in matches[matchpthr][query_id]:
                    matches[matchpthr][query_id]['hmmstart'] = []
                if 'hmmend' not in matches[matchpthr][query_id]:
                    matches[matchpthr][query_id]['hmmend'] = []

                matches[matchpthr][query_id]['hmmstart'].append(mark[6])
                matches[matchpthr][query_id]['hmmend'].append(mark[7])

                # print(matches)

            elif m.match('\s+==') and align_found_n == 1:
                # print("INSIDE IF 5: {}".format(line.strip()))
                # query_id = m.group(1)
                # print(query_id)
                # print(align_found_n)

                line = fp.readline()

                if 'hmmalign' not in matches[matchpthr][query_id]:
                    matches[matchpthr][query_id]['hmmalign'] = []

                matches[matchpthr][query_id]['hmmalign'].append(line.split()[2])

                line = fp.readline()

                line = fp.readline()
                if 'matchalign' not in matches[matchpthr][query_id]:
                    matches[matchpthr][query_id]['matchalign'] = []

                matches[matchpthr][query_id]['matchalign'].append(
                    line.split()[2])

                line = fp.readline()

            elif m.match('\A>> (\S+)'):
                align_found_n += 1
                # print("INSIDE IF 3: {}".format(line.strip()))
                matchpthr = m.group(1)
                # print(matchpthr)
                # print(align_found_n)


            elif m.match('\A\/\/'):
                # print("END BLOCK")
                align_found_n = 0
                # break
            # else:
            #     print("NOT MATCHED: {}".format(line.strip()))

            line = fp.readline()

    fp.close()
    return(matches)


def parsehmmsearch(hmmer_out):

    matches = {}

    with open(hmmer_out) as fp:
        score_store = {}

        match_store = {}

        # match_store['POPTR|EnsemblGenome=POPTR_0018s04850|UniProtKB=B9INH6'] = {}
        # match_store['POPTR|EnsemblGenome=POPTR_0018s04850|UniProtKB=B9INH6']['score'] = 900
        # match_store['POPTR|EnsemblGenome=POPTR_0018s04850|UniProtKB=B9INH6']['panther_id'] = 'TEST_PNTR_ID'

        store_align = 0

        matchpthr = None
        query_id = None

        line = fp.readline()
        while line:
            m = re_matcher(line)
            # print("Line: {}".format(line.strip()))

            if line.startswith('#') or not line.strip():
                # print("INSIDE IF 1: {}".format(line.strip()))
                line = fp.readline()
                continue
            elif m.match('\AQuery:\s+(PTHR[0-9]+)'):
                # print("INSIDE IF 2: {}".format(line.strip()))
                matchpthr = m.group(1)
                # print(matchpthr)

                fp.readline()
                fp.readline()
                fp.readline()
                fp.readline()
                line = fp.readline()
                # print(line)

                while line.strip():
                    m = re_matcher(line)
                    if m.match('\s+------\sinclusion\sthreshold'):
                        # print("INCLUSION THRESHOLD:")
                        # print(line)
                        break
                    
                    # print("STRIP:")
                    # print(line)
                    
                    score_array = line.split()

                    # curr_query_id = score_array[8]
                    # # print(curr_query_id)
                    # curr_score = float(score_array[1])
                    # # print(curr_score)

                    score_store[score_array[8]] = float(score_array[1])


                    line = fp.readline()

                # print("END IF 2:")
                # print(line)

            elif m.match('\A>>\s+(\S+)'):
                # print("INSIDE IF 3: {}".format(line.strip()))
                query_id = m.group(1)
                store_align = 0

                if not query_id in score_store:
                    # print("QUERY ID UNDER THRESHOLD")
                    line = fp.readline()
                    continue
    
                # print(query_id)
                if query_id not in match_store or score_store[query_id] > match_store[query_id]['score']:
                    # print("NEED TO STORE MATCH")
                    store_align = 1
                    # if query_id not in match_store:
                    match_store[query_id] = {
                            'panther_id': matchpthr, 'score': score_store[query_id], 'align': {
                                'hmmstart': [], 'hmmend': [], 'hmmalign': [], 'matchalign': []
                            } }



            elif m.match('\s+==') and store_align:
                # print("INSIDE IF 4: {}".format(line.strip()))

                line = fp.readline()
                hmmalign_array = line.split()

                match_store[query_id]['align']['hmmstart'].append(hmmalign_array[1])
                match_store[query_id]['align']['hmmend'].append(hmmalign_array[3])
                match_store[query_id]['align']['hmmalign'].append(hmmalign_array[2])

                line = fp.readline()

                line = fp.readline()

                match_store[query_id]['align']['matchalign'].append(line.split()[2])

                line = fp.readline()

            elif m.match('\A\/\/'):
                # print("END BLOCK")
                # print(match_store)
                score_store = {}
                # BREAK FOR DEBUG
                # break

            # else:
            #     print("NOT MATCHED: {}".format(line.strip()))

            line = fp.readline()

    fp.close()

    for query_id in match_store:
        panther_id = match_store[query_id]['panther_id']

        if panther_id not in matches:
            matches[panther_id] = {}

        matches[panther_id][query_id] = match_store[query_id]['align']


    # print(matches)
    return(matches)




def get_args():
    """
    Command line arguments parser.
    """

    ap = argparse.ArgumentParser(
        prog='treegrafter.py', description="TreeGrafter",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ap.add_argument('-f', '--fasta', required=True,
        help="input fasta file")

    ap.add_argument(
        '-hb', '--hbin', default=None,
        help='path to hmmr bin (default PATH)')

    ap.add_argument(
        '-hm', '--hmode', default='hmmscan', choices=['hmmscan', 'hmmsearch'],
        help='hmmr mode to use')

    ap.add_argument(
        '-hc', '--hcpus', default=1, type=int,
        help="number of hmmr cpus")

    ap.add_argument(
        '-ho', '--hout', default=None,
        help="existing hmmr output file")

    ap.add_argument(
        '-o', '--out', required=True,
        help='output file name')

    ap.add_argument(
        '-d', '--data', required=True,
        help='panther data directory')

    ap.add_argument(
        '-t', '--tmp',
        help='tmp work directory (default tmp folder in provided panther data directory)')




    # ap.add_argument(
    #     '-rb', '--rbin', default='',
    #     help='path to RAxML bin (default PATH)')


    # parser.add_argument('-l', '--log', type=str, default=None, help='log file')
    # parser.add_argument("-ll", "--logLevel", default="ERROR",
    #                     choices=["DEBUG", "INFO",
    #                              "WARNING", "ERROR", "CRITICAL"],
    #                     help='log level filter. All levels <= choice will be displayed')

    args = vars(ap.parse_args())


    # print(args)
    return args


def align_lenght(pthr):
    pthr_fasta_file = options['msf_tree_folder'] + pthr + '.AN.fasta'

    try:
        with open(pthr_fasta_file) as f:
            file = f.read().split('>')
            first_seq = file[1]
            first_seq = re.sub('\A[^\n]+', '', first_seq)
            first_seq = re.sub('\n', '', first_seq)
            seq_lenght = len(first_seq)
    except IOError:
        print("Could not find alignment for " + pthr)
        return(0)

    return(seq_lenght)


def get_annotations():
    annot_file = options['data_folder'] + 'PAINT_Annotations/PAINT_Annotatations_TOTAL.txt'
    # print(annot_file)

    annot = {}
    with open(annot_file) as f:
        for line in f:
            line = line.strip()
            # print(line)
            line_array = line.split("\t")
            # print(line_array)
            annot[line_array[0]] = line_array[1]

    return annot





if __name__ == '__main__':

    args = get_args()

    global options
    options = {}
    options['data_folder'] = os.path.abspath(args['data']) + '/'
    options['fasta_input'] = args['fasta']
    options['out_file'] = args['out']
    options['hmmr_mode'] = args['hmode']
    options['hmmr_bin'] = args['hbin']
    options['hmmr_cpus'] = args['hcpus']
    options['hmmr_out'] = args['hout']

    options['msf_tree_folder'] = options['data_folder'] + 'Tree_MSF/'
    if args['tmp'] is None:
        options['tmp_folder'] = options['data_folder'] + 'tmp/'
    else:
        options['tmp_folder'] = args['tmp'] + '/'
    
    # print(json.dumps(options, indent=4))
    

    # testing
    # options['msf_tree_folder'] = '/home/tgrego/dev/treegrafter/Test/PANTHER_mini/PANTHER12_Tree_MSF/'

    # testing

    
    # hmmsearch_file = '/home/tgrego/dev/treegrafter/Test/sample.fasta.hmmsearch.out'
    # hmmscan_file = '/home/tgrego/dev/treegrafter/Test/sample.fasta.hmmscan.out'

    # matches = parsehmmsearch(hmmsearch_file)
    # matches = parsehmmscan(hmmscan_file)
    # print(matches)

    if options['hmmr_out'] is None:
        runhmmr()

    matches = parsehmmr(options['hmmr_out'])

    annotations = get_annotations()
    # print(annotations)

    results = []

    for pthr in matches:
        print(pthr)
        logging.info('Processing panther id ' + pthr + "...\n")

        pthr_align_lenght = align_lenght(pthr)


        for query_id in matches[pthr]:
            # print(query_id)
            query_id_str = stringify(query_id)
            # print(query_id_str)
            query_msf = _querymsf(matches[pthr][query_id], pthr_align_lenght)
            # print(query_msf)

            fasta_file = _generateFasta(pthr, query_id_str, query_msf)
            # print(fasta_file)

            result_string = _run_raxml(pthr, query_id_str, fasta_file, annotations)
            results.append(result_string)
            # print(result_string)

    # print(results)
    file = open(options['out_file'], 'w')
    for line in results:
        file.write(line)
    file.close()







    # test_input_querymsf = [
    #     {
    #         'hmmalign': [
    #             'klialDlDGTLlnskkeiskrtlealkeakerGvkvviaTGrsraaviellkeldlgsplvtlnGalvyskqgevlfernldpevlrellelaeeegvalvaysedrssplveslhtiykepkvekvesleklleeapiskvlflstdeeklealrevleealegelsvtrsapdfleivpkgvsKgsglkrlleelgisleeviafGDgeNDlemLelaglgvamgnasekvkevadvvtasndedGvakaleky'
    #         ],
    #         'hmmend': [
    #             '263'
    #         ],
    #         'hmmstart': [
    #             '8'
    #         ],
    #         'matchalign': [
    #             'LVIFTDIDGTLY-GDFH----IHEAFKRFITNGLFLVYSTGRNLQSFKDLQKNVHLPDILVGSCGSEIYQLGQDEFETNPYNQN--QAWIQYITQDNWDLQALYD-----------FVKKEFP----AAWPNLSEGVSLYKGSFLLTDSRQRDKLDVLMKKAFLNKYIISGHGHRFLDILPERADKGSSLQFVCKILKTDYTKSAAFGDSLNDVDLLCCAGQGFIVANAQSQNQRFQNVKVSYHEGDAIAKYLQQI'
    #         ]
    #     }, {
    #         'hmmalign': [
    #             'iklialDlDGTLlnskkeiskrtlealkeakerGvkvviaTGrsraaviellkeldlgsplvtlnGalvyskqgevlfernldpevlrellelaeeegvalvaysedrssplveslhtiykepkvekvesleklleeapiskvlflstdeeklealrevleealegelsvtrsapdfleivpkgvsKgsglkrlleelgisleeviafGDgeNDlemLelaglgvamgnasekvkevadvvtasndedGvakaleky'
    #         ],
    #         'hmmend': [
    #             '263'
    #         ],
    #         'hmmstart': [
    #             '7'
    #         ],
    #         'matchalign': [
    #             'YRVFVFDLDGTLLNDNLEISEKDRRNIEKL-SRKCYVVFASGRMLVSTLNVEKKFKRTFPTIAYNGAIVYLPEEGVILNEKIPPEVAKDIIEYIKPLNVHWQAYIDDV---LSKDNEKSYARHSYRVEPNLSELVSKMGTTKLLLIDT-PERLDELKEILSERFKDVVKVFKSFPTYLEIVPKNVDKGKALRFLRERMNWKKEEIVVFGDNENDLFMFEEAGLRVAMENAIEKVKEASDIVTLTNNDSGVSYVLERI'
    #         ]
    #     }
    #     # , {
    #     #     'hmmalign': [
    #     #         'kpkiklialDlDGTLlnskkeiskrtlealkeakerGvkvviaTGrsraaviellkeldlgsplvtlnGalvyskqgevlfernldpevlrellelaeeegvalvaysedrssplveslhtiykepkvekvesleklleeapiskvlflstdeeklealrevleealegelsvtrsapdfleivpkgvsKgsglkrlleelgisleeviafGDgeNDlemLelaglgvamgnasekvkevadvvtasndedGvakalekyll',
    #     #         'kpkik'
    #     #     ],
    #     #     'hmmend': [
    #     #         '265',
    #     #         '275'
    #     #     ],
    #     #     'hmmstart': [
    #     #         '4',
    #     #         '270'
    #     #     ],
    #     #     'matchalign': [
    #     #         'RMPIKAVVTDLDGTLLDPQHCISNYAAEVLKKIKEKGICFIVATGRPYAEVFNRIRHCHLPDYIITSNGARIHDGAFNVVREHNLRPELVESLARVRTVKDPATNIYRGLT--PEVSAFHTDFQCTDRERFYELQ-ASELGDVHEIWFAGD-HDELVLLDNALREKYPGDLCCTFSLPHLLDCVPAGVNKGNGVREAAEMLGLALDEVACFGDGMNDESMLQVTSTSFIMANAQQRLKAVPHAQIISNADDGVAKKLEEMFF',
    #     #         'RMPIK'
    #     #     ]
    #     # }
    # ]

    # test_query_id_list = iter([
    #     'TETTS|EnsemblGenome=TTHERM_00590270|UniProtKB=I7M2B5',
    #     'THEMA|EnsemblGenome=TM_0651|UniProtKB=Q9WZB9',
    #     'LEIMA|EnsemblGenome=LmjF.28.1370|UniProtKB=Q4Q8C9',
    #     'BACTN|EnsemblGenome=BT_4131|UniProtKB=Q8A090',
    #     'SYNY3|Gene=BAA18460|UniProtKB=P74365'
    # ])

    # test_lenght = iter([
    #     266,
    #     266,
    #     280,
    #     266,
    #     266
    # ])

    # for test_data in test_input_querymsf:
    #     test_query_id = next(test_query_id_list)

    #     test_query_id_str = stringify(test_query_id)


    #     query_msf = _querymsf(test_data, next(test_lenght))

    #     print(query_msf)

    #     fasta_file = _generateFasta('PTHR10000', test_query_id_str, query_msf)

    #     print(fasta_file)

    #     result_string = _run_raxml('PTHR10000', test_query_id_str, fasta_file)

    #     print(result_string)



    
