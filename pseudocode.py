

def graftMatches(matches):
    # matches is a dictionary containing as
    # keys panther family matches (eg 'PTHR10000') and as
    # values a dictionary with the matched seq query_id (eg 
    # 'ARATH|TAIR=AT2G25870|UniProtKB=Q8L5Z4') and values
    # the hmmer output (hmmstart, hmmend, hmmalign and matchalign, all as arrays) 

    for pathr in matches.keys():

        next if annotationsNotAvailable

        pthrAlignLength = get_alignment_lenght  # from file "pantherdir/pathr_id.AN.fasta"

        for query_id, matchdata in matches.pathr.items():
            # query_id=>matchdata is a match sub dictinary for a panther family
            querymsf = _querymsf(query_id, matchdata, pthrAlignLength)

            # queryfasta = _generateFasta()



def _querymsf(query_id, matchdata, pthrAlignLength):
    # matchdata contains: hmmstart, hmmend, hmmalign and matchalign, as arrays (multiple values possible)

    querymsf = '' # placeholder for the return sequence

    # N-terminaly padd the sequence
    # position 1 until start is filled with '-'
    start = matchdata[hmmstart][0]
    querymsf .= (start - 1 ) x '-'


    #For the first element/domain, extract the query string
    matchalign = matchdata[matchalign][0]
    hmmalign = matchdata[hmmalign][0]
    for i in range(0, len(hmmalign)):
        next if (hmmalign[i:i+1] eq "."))  # hmm insert state
        querymsf .=matchalign[i:i+1]


    #Now we need to add in any additional domains in the next part
    # this would be where multiple values occur
    for j in range(1, len(matchdata[matchalign])):
        start = matchdata[hmmstart][j]
        end = matchdata[hmmend][j-1]
        #This bridges the gap between the hits
        querymsf .= (start - end) x '-'

        # repeat logic from first element/domain
        # refactor as a function to reuse
        matchalign=matchdata[matchalign][j]
        hmmalign=matchdata[hmmalign][j]
        for i in range(0, len(hmmalign)):
            next if (hmmalign[i:i+1] eq "."))  # insert state, so skipping
            querymsf .= matchalign[i:i+1]


    # get the end of the last element/domain
    prevend=matchdata[hmmend][-1]
    # and padd out to fill the msf lenght
    querymsf .= (pthrAlignLength -1 - prevend) x '-'

    if (len(querymsf) != pthrAlignLength) :
        # then something is wrong
        return 0
    
    querymsf=querymsf.upper()

    return querymsf


# def _generateFasta():





def main:

    args = process_args()

    hmmer_results = runhmmer() # or use stored result file

    matches = parsehmmer(hmmer_results)

    results = graftMatches(matches)

    print_results(results)
