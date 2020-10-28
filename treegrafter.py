#!/usr/bin/env python

import re



def _querymsf(query_id, matchdata, pthrAlignLength):
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
            # This bridges the gap between the hits
            querymsf += (start - end) * '-'

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
        return 0

    return querymsf.upper()


def _generateFasta(pathr, query_id, querymsf):

    # stringify query_id
    query_id = re.sub('[^\w]', '_', query_id)

    # use static dir paths for testing. these are provided
    fasta_in_dir = '/home/tgrego/dev/treegrafter/Test/PANTHER_mini/PANTHER12_Tree_MSF/'
    fasta_out_dir = '/home/tgrego/dev/treegrafter/Test/PANTHER_mini/tmp/'
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



if __name__ == '__main__':

    test_input_querymsf = [
        {
            'hmmalign': [
                'klialDlDGTLlnskkeiskrtlealkeakerGvkvviaTGrsraaviellkeldlgsplvtlnGalvyskqgevlfernldpevlrellelaeeegvalvaysedrssplveslhtiykepkvekvesleklleeapiskvlflstdeeklealrevleealegelsvtrsapdfleivpkgvsKgsglkrlleelgisleeviafGDgeNDlemLelaglgvamgnasekvkevadvvtasndedGvakaleky'
            ],
            'hmmend': [
                '263'
            ],
            'hmmstart': [
                '8'
            ],
            'matchalign': [
                'LVIFTDIDGTLY-GDFH----IHEAFKRFITNGLFLVYSTGRNLQSFKDLQKNVHLPDILVGSCGSEIYQLGQDEFETNPYNQN--QAWIQYITQDNWDLQALYD-----------FVKKEFP----AAWPNLSEGVSLYKGSFLLTDSRQRDKLDVLMKKAFLNKYIISGHGHRFLDILPERADKGSSLQFVCKILKTDYTKSAAFGDSLNDVDLLCCAGQGFIVANAQSQNQRFQNVKVSYHEGDAIAKYLQQI'
            ]
        }, {
            'hmmalign': [
                'iklialDlDGTLlnskkeiskrtlealkeakerGvkvviaTGrsraaviellkeldlgsplvtlnGalvyskqgevlfernldpevlrellelaeeegvalvaysedrssplveslhtiykepkvekvesleklleeapiskvlflstdeeklealrevleealegelsvtrsapdfleivpkgvsKgsglkrlleelgisleeviafGDgeNDlemLelaglgvamgnasekvkevadvvtasndedGvakaleky'
            ],
            'hmmend': [
                '263'
            ],
            'hmmstart': [
                '7'
            ],
            'matchalign': [
                'YRVFVFDLDGTLLNDNLEISEKDRRNIEKL-SRKCYVVFASGRMLVSTLNVEKKFKRTFPTIAYNGAIVYLPEEGVILNEKIPPEVAKDIIEYIKPLNVHWQAYIDDV---LSKDNEKSYARHSYRVEPNLSELVSKMGTTKLLLIDT-PERLDELKEILSERFKDVVKVFKSFPTYLEIVPKNVDKGKALRFLRERMNWKKEEIVVFGDNENDLFMFEEAGLRVAMENAIEKVKEASDIVTLTNNDSGVSYVLERI'
            ]
        }, {
            'hmmalign': [
                'kpkiklialDlDGTLlnskkeiskrtlealkeakerGvkvviaTGrsraaviellkeldlgsplvtlnGalvyskqgevlfernldpevlrellelaeeegvalvaysedrssplveslhtiykepkvekvesleklleeapiskvlflstdeeklealrevleealegelsvtrsapdfleivpkgvsKgsglkrlleelgisleeviafGDgeNDlemLelaglgvamgnasekvkevadvvtasndedGvakalekyll',
                'kpkik'
            ],
            'hmmend': [
                '265',
                '275'
            ],
            'hmmstart': [
                '4',
                '270'
            ],
            'matchalign': [
                'RMPIKAVVTDLDGTLLDPQHCISNYAAEVLKKIKEKGICFIVATGRPYAEVFNRIRHCHLPDYIITSNGARIHDGAFNVVREHNLRPELVESLARVRTVKDPATNIYRGLT--PEVSAFHTDFQCTDRERFYELQ-ASELGDVHEIWFAGD-HDELVLLDNALREKYPGDLCCTFSLPHLLDCVPAGVNKGNGVREAAEMLGLALDEVACFGDGMNDESMLQVTSTSFIMANAQQRLKAVPHAQIISNADDGVAKKLEEMFF',
                'RMPIK'
            ]
        }
    ]

    test_query_id_list = iter([
        'TETTS|EnsemblGenome=TTHERM_00590270|UniProtKB=I7M2B5',
        'THEMA|EnsemblGenome=TM_0651|UniProtKB=Q9WZB9',
        'LEIMA|EnsemblGenome=LmjF.28.1370|UniProtKB=Q4Q8C9',
        'BACTN|EnsemblGenome=BT_4131|UniProtKB=Q8A090',
        'SYNY3|Gene=BAA18460|UniProtKB=P74365'
    ])

    test_lenght = iter([
        266,
        266,
        280,
        266,
        266
    ])


    for test_data in test_input_querymsf:
        test_query_id = next(test_query_id_list)


        query_msf = _querymsf(test_query_id, test_data, next(test_lenght))

        # print(query_msf)

        fasta_file = _generateFasta('PTHR10000', test_query_id, query_msf)
