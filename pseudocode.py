

def graftMatches(matches):
    # matches is a dictionary containing as
    # keys panther family matches (eg 'PTHR10000') and as
    # values a dictionary with the matched seq query_id (eg 
    # 'ARATH|TAIR=AT2G25870|UniProtKB=Q8L5Z4') and values
    # the hmmer output (hmmstart, hmmend, hmmalign and matchalign, all as arrays) 
    
    allResults = []

    for pathr in matches.keys():

        next if annotationsNotAvailable

        pthrAlignLength = get_alignment_lenght  # from file "pantherdir/{pathr}.AN.fasta"

        for query_id, matchdata in matches.pathr.items():
            # query_id=>matchdata is a match sub dictinary for a panther family
            querymsf = _querymsf(query_id, matchdata, pthrAlignLength)

            queryfasta = _generateFasta(pathr, query_id, querymsf)

            result_string = _runRAxMLAndAnnotate(pathr, queryid, queryfasta)

            allResults.append(result_string) if result_string

    return allResults


def _querymsf(query_id, matchdata, pthrAlignLength):
    # matchdata contains: hmmstart, hmmend, hmmalign and matchalign, as arrays (multiple values possible)

    querymsf = '' # placeholder for the return sequence

    # N-terminaly padd the sequence
    # position 1 until start is filled with '-'
    start = matchdata[hmmstart][0]
    querymsf .= (start - 1 ) * '-'


    #For the first element/domain, extract the query string
    matchalign = matchdata[matchalign][0]
    hmmalign = matchdata[hmmalign][0]
    for i in range(0, len(hmmalign)):
        next if (hmmalign[i:i+1] eq "."))  # hmm insert state
        querymsf .= matchalign[i:i+1]


    #Now we need to add in any additional domains in the next part
    # this would be where multiple values occur
    for j in range(1, len(matchdata[matchalign])):
        start = matchdata[hmmstart][j]
        end = matchdata[hmmend][j-1]
        #This bridges the gap between the hits
        querymsf .= (start - end) * '-'

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
    querymsf .= (pthrAlignLength -1 - prevend) * '-'

    if (len(querymsf) != pthrAlignLength) :
        # then something is wrong
        return 0
    
    querymsf=querymsf.upper()

    return querymsf


def _generateFasta(pathr, query_id, querymsf):

    # replace special chars with underscore
    query_id=re.sub('[^\w]', '_', query_id)

    # write out query fasta
    queryfasta = open(args[directory]."/tmp/{queryid}.{pathr}.fasta", 'w')
    queryfasta.write('>' . query_id . "\n")
    queryfasta.write(querymsf . "\n") # original split into 80 char lines

    # read in align msf
    hmmalignmsf = open(args[pantherdir]."/{pathr}.AN.fasta", 'r')
    # and write out to query fasta
    for line in hmmalignmsf.readlines():
        queryfasta.write(line)

    return queryfasta


def _runRAxMLAndAnnotate(pathr, queryid, queryfasta):

    # replace special chars with underscore
    # same as in _generateFasta(), factor out code duplication
    query_id=re.sub('[^\w]', '_', query_id)

    bifurnewick=args[pantherdir]."/{panthr}.bifurcate.newick"
    if no bifurnewick dir exists:
        print("no bifurcate newickfile for {pathr}}\n")
        return 0

    raxmldir=args[directory]."/tmp/{panthr}_{queryid}_raxml{PID:process id}"
    os.mkdir(raxmldir)

    # raxml-ng using epa-ng is under development
    raxml_cmd=args[raxmlloc]."raxmlHPC-PTHREADS-SSE3 -f y -p 12345 -t {bifurnewick} -G 0.05 -m PROTGAMMAWAG -T 4 -s {queryfasta} -n $matchpthr -w {raxmldir} >/dev/null"
    try:
        os.system(raxml_cmd)
    except Exception as e:
        print("caught error for {pathr}, {queryid}:" e)

    ## now find the location of the graft sequence in the tree
    mapANs=mapto(raxmldir, pathr, queryid)

    if mapANs == 'root':
        annotation = args[annotations][{pathr}.":".{mapANs}]
    else:
        commonAN = commonancestor(mapANs, pathr)
        annotation = args[annotations][{pathr}.":".{commonAN}]

    result="{queryid}\t{pathr}\t{annotation}\n"

    return result




# this function retrieves the location of the query sequence in the raxml output
def mapto(raxmldir, pathr, queryid):
...
returns childrenids if any or 'root'


# find common ancestor of a list of ANs which are leaf genes
def commonancestor(mapANs, pathr)
...
return commonancestor if any or 'root'


def main:

    args = process_args()

    hmmer_results = runhmmer() # or use stored result file

    matches = parsehmmer(hmmer_results)

    results = graftMatches(matches)

    print_results(results)
