import Levenshtein
from spoa import poa
import time
from multiprocessing import Pool
import random
import argparse
import math

def NW(s,anchor):
    row = len(s)+1
    col = len(anchor)+1
    matrix = [[0]*col for i in range(row)]
    for j in range(col):
        matrix[0][j] = -j
    # print(matrix)
    for i in range(1,row):
        for j in range(1,col):
            if s[i-1] == anchor[j-1]:
                match = matrix[i-1][j-1] + 1
            else:
                match = matrix[i-1][j-1] - 1
            delete = matrix[i-1][j] - 1
            insert = matrix[i][j-1] - 1
            matrix[i][j] = max(match,delete,insert)
    # for i in matrix:
    #     print(i)
    max_score = -1
    end = -1
    for i in range(row):
        if matrix[i][-1] > max_score:
            max_score = matrix[i][-1]
            end = i
    # print(end)
    alignment_s = []
    alignment_a = []
    i = end
    j = len(anchor)
    while i > 0 and j > 0:
        if (i > 0) and (j > 0) and matrix[i][j] == matrix[i-1][j-1] + 2*(s[i-1] == anchor[j-1]) - 1:
            alignment_s.append(s[i-1])
            alignment_a.append(anchor[j-1])
            i -= 1
            j -= 1
        elif (i > 0) and matrix[i][j] == matrix[i-1][j] - 1:
            alignment_s.append(s[i-1])
            alignment_a.append('-')
            i -= 1
        else:
            alignment_s.append('-')
            alignment_a.append(anchor[j-1])
            j -= 1
    # for x in matrix:
    #     print(x)
    return s[:i]+(''.join(reversed(alignment_s))) + s[end:],'-'*i+''.join(reversed(alignment_a))+(len(s)-end)*'-'#, i, end-1
def ORIGINAL_DoubleSidedGreedyMedian(cluster):
    overlap = 2
    strand_length = STRAND_LENGTHS[0]
    rev = list(map(lambda x: x[len(x)-1:len(x)//2-overlap:-1], cluster))
    cluster = list(map(lambda x: x[:len(x)//2+1+overlap],cluster))
    LtoR = Levenshtein.median(cluster)
    RtoL = Levenshtein.median(rev)
    if len(LtoR) + len(RtoL) < strand_length:
        diff = strand_length - len(LtoR) - len(RtoL)
        return LtoR + LtoR[-1]*diff + RtoL
    return LtoR + RtoL[(strand_length-len(LtoR))-1::-1]
def DoubleSidedGreedyMedian(cluster):
    rev = list(map(lambda x: x[len(x)-1:len(x)//2-overlap:-1], cluster))
    cluster = list(map(lambda x: x[:len(x)//2+1+overlap],cluster))
    LtoR = Levenshtein.median(cluster)
    RtoL = Levenshtein.median(rev)
    if len(LtoR) + len(RtoL) < strand_length:
        diff = strand_length - len(LtoR) - len(RtoL)
        return LtoR + LtoR[-1]*diff + RtoL
    return LtoR + RtoL[(strand_length-len(LtoR))-1::-1]
def splitByAnchor(s, anchor):
    n = len(s)
    middle = s[(n//2)-(D//2):(n//2)-(D//2)+D]
    alignment = NW(middle, anchor)
    # print(s, anchor)
    # print(((len(s)//2)-(D//2))*' '+alignment[0])
    # print(((len(s)//2)-(D//2))*' '+alignment[1])

    # alignment = ('GCAG-AGTCTGGTAG', '-CAGTAG-CT-----')
    if alignment[1] == '-'*(len(middle)+1):
        return ('','')
    if alignment[1][0] != '-':
        c = 0
        while alignment[0][c] == '-':
            c += 1
        l = s[:(n//2)-(D//2)-c]
    else:
        c = 0
        while alignment[1][c] == '-':
            c += 1
        l = s[:(n//2)-(D//2)+c]
    if alignment[1][-1] != '-':
        c = 1
        while alignment[0][-c] == '-':
            c += 1
        c -= 1
        r = s[(n//2)-(D//2)+D+c:]
    else:
        c = 1
        while alignment[1][-c] == '-':
            c+=1
        c -= 1
        r = s[(n//2)-(D//2)+D-c:]
    # print(s)
    # print(' '*((n//2)-(D//2)) + middle +' '*(55-(D//2)))
    # print(' '*(55-(D//2)) + alignment[0])
    # print(' '*(55-(D//2)) + alignment[1])
    # print(l)
    # print(' '*(n-len(r))+r)
    return l, r
def AnchorReconstruct(cluster):
    strand_length = STRAND_LENGTHS[1]
    L, R = [], []
    for i in range(len(cluster)):
        # print(i)
        l, r = splitByAnchor(cluster[i], anchor)
        if l and r:
            L.append(l)
            R.append(r)
    # for i in range(len(cluster)):
    #     print(cluster[i])
    #     print(L[i])
    #     print(' '*(len(cluster[i])-len(R[i]))+R[i])
    #     print()
    return f(L) +  f(R)
bases = ['A', 'C', 'G', 'T']
blank = '-'
window_size = 3
overlap = 2
def majority(list):
    chars = set(list)
    chars.discard(blank)
    return max(chars, key=list.count)
def refine_majority(clu, i):
    ems = []
    for ele in clu:
        if len(ele) > i:
            ems.append(ele[i])
    if len(ems) == 0:
        ems.append(random.choice(bases))
    return ems
def NW_recover_strand(cluster, strand_len):
    c, alignment = poa(cluster, algorithm=1)
    ref_size = strand_len
    length = len(alignment[0])
    pos_to_unkown_num = {}
    for pos in range(length):
        unkown_num = 0
        for strand in alignment:
            if strand[pos] == '-':
                unkown_num = unkown_num + 1
        pos_to_unkown_num[pos] = unkown_num
    a = sorted(pos_to_unkown_num.items(), key=lambda kv:(kv[1]), reverse = True)
    skip_index_num = length - ref_size
    skip_index = []
    loop = 0
    for pos, num in a:
        if loop == skip_index_num:
            break
        skip_index.append(pos)
        loop = loop + 1
    new_cluster = []
    for i in range(len(alignment)):
        new_cluster.append("")
    for pos in range(0, length):
        if pos in skip_index:
            continue
        for id in range(len(alignment)):
            new_cluster[id] = new_cluster[id]+alignment[id][pos]
    
    ans = ''
    recovered = ''
    for i in range(0, strand_len - 1):
        ch = majority(refine_majority(new_cluster, i))
        recovered += ch
        ans = ans[:0] + recovered
    last_ch = majority(refine_majority(new_cluster, strand_len - 1))
    ans += last_ch

    return ans
def recover_strand(cluster, strand_len):
    ans = ''
    recovered = ''

    for i in range(0, strand_len - 1):
        ch = majority(refine_majority(cluster, i))

        for j in range(len(cluster)):

            if len(cluster[j]) == i:
                cluster[j] += ch

            if cluster[j][i] != ch:

                ch2 = majority(refine_majority(cluster, i + 1))

                ch3_flag = -1
                if i + 2 < strand_len:
                    ch3_flag = 1
                    ch3 = majority(refine_majority(cluster, i + 2))

                ch4_flag = -1
                if i + 3 < strand_len:
                    ch4_flag = 1
                    ch4 = majority(refine_majority(cluster, i + 3))

                ch5_flag = -1
                if i + 4 < strand_len:
                    ch5_flag = 1
                    ch5 = majority(refine_majority(cluster, i + 4))

                if len(cluster[j]) > i + 2:
                    if cluster[j][i] == ch2 and (ch3_flag == -1 or cluster[j][i + 1] == ch3):  # erasure error
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i:]

                    elif cluster[j][i + 1] == ch and cluster[j][i + 2] == ch2:  # insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i + 1:]

                    elif cluster[j][i + 1] == ch2 and (ch3_flag == -1 or cluster[j][i + 2] == ch3):  # subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i + 1:]

                    elif cluster[j][i + 1] != ch2:

                        if cluster[j][i] == ch3 and (ch4_flag == -1 or cluster[j][i + 1] == ch4):  # erasure
                            cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i:]

                        elif len(cluster[j]) > i + 3:
                            if cluster[j][i + 2] == ch3 and (ch4_flag == -1 or cluster[j][i + 3] == ch4):  # subs
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i + 1:]

                            elif cluster[j][i + 2] == ch and cluster[j][i + 3] == ch2:  # insertion
                                cluster[j] = cluster[j][:i] + cluster[j][i + 2:]

                            elif cluster[j][i + 1] == ch3 and (ch4_flag == -1 or cluster[j][i + 2] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i + 1:]

                            elif cluster[j][i + 1] == ch3 and (ch4_flag == -1 or cluster[j][i + 3] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i + 1:]

                            elif cluster[j][i + 2] == ch2 and cluster[j][i + 3] == ch3:
                                cluster[j] = cluster[j][:i] + ch + cluster[j][i + 2:]

                            elif cluster[j][i] == ch3 and (ch4_flag == -1 or cluster[j][i + 3] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + ch3 + cluster[j][i + 3:]

                            elif cluster[j][i + 2] != ch3:

                                if cluster[j][i] == ch4 and (ch5_flag == -1 or cluster[j][i + 1] == ch5):  # erasure
                                    cluster[j] = cluster[j][:i] + ch + ch2 + ch3 + cluster[j][i:]

                                elif len(cluster[j]) > i + 4:
                                    if cluster[j][i + 3] == ch4 and (
                                            ch5_flag == -1 or cluster[j][i + 4] == ch5):  # subs
                                        cluster[j] = cluster[j][:i] + ch + ch2 + ch3 + cluster[j][i + 1:]

                                    elif cluster[j][i + 3] == ch and cluster[j][i + 4] == ch2:  # insertion
                                        cluster[j] = cluster[j][:i] + cluster[j][i + 3:]

                elif len(cluster[j]) == i + 2:
                    if cluster[j][i] == ch2:  # erasure error
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i:]

                    elif cluster[j][i + 1] == ch2:  # subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i + 1:]

                    elif cluster[j][i + 1] == ch:  # insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i + 1:]

                    else:
                        cluster[j] = cluster[j][:i] + ch

        recovered += ch
        ans = ans[:0] + recovered

    last_ch = majority(refine_majority(cluster, strand_len - 1))
    ans += last_ch

    return ans
def clean_up(file_params, strand_num, myPath):
    if myPath and os.path.exists(myPath):
       shutil.rmtree(myPath)

    # if cleanup_error==1:
    #     if os.path.exists("results/error" + file_params + ".txt"):
    #         os.remove("results/error" + file_params + ".txt")
    #     if os.path.exists("results/error" + file_params + "_dp" + ".txt"):
    #         os.remove("results/error" + file_params + "_dp" + ".txt")
def ANCHORED_reconstruct(cluster):
    strand_length = STRAND_LENGTHS[1]
    strand_num=len(cluster)
    
    file_params = "w" + str(window_size) + "n" + str(strand_num) + "l" + str(strand_length)

    rev_cluster = []
    myPath=""

    for i in range(0, len(cluster)):
        rev_cluster.append(cluster[i][::-1])

    # Two sided BMA
    if ALG==0:
        mj = recover_strand(cluster, strand_length)
        rev_mj = recover_strand(rev_cluster, strand_length)
        rev_rev_mj = rev_mj[::-1]
        mj = mj[0:int(strand_length / 2) - 1] + rev_rev_mj[int(strand_length / 2) - 1:strand_length]


    # Single sided BMA
    elif ALG==1:
        mj = recover_strand(cluster, strand_length)
    # Needleman-Wunsch
    else:
        mj = NW_recover_strand(cluster, strand_length)
    #clean_up(file_params, strand_num, myPath)
    return mj
def ORIGINAL_reconstruct(cluster):
    strand_length = STRAND_LENGTHS[0]
    strand_num=len(cluster)
    
    file_params = "w" + str(window_size) + "n" + str(strand_num) + "l" + str(strand_length)

    rev_cluster = []
    myPath=""

    for i in range(0, len(cluster)):
        rev_cluster.append(cluster[i][::-1])

    # Two sided BMA
    if ALG==0:
        mj = recover_strand(cluster, strand_length)
        rev_mj = recover_strand(rev_cluster, strand_length)
        rev_rev_mj = rev_mj[::-1]
        mj = mj[0:int(strand_length / 2) - 1] + rev_rev_mj[int(strand_length / 2) - 1:strand_length]


    # Single sided BMA
    elif ALG==1:
        mj = recover_strand(cluster, strand_length)
    # Needleman-Wunsch
    else:
        mj = NW_recover_strand(cluster, strand_length)
    #clean_up(file_params, strand_num, myPath)
    return mj
def reconstruct_clusters(clusters):
   pool = mp.Pool(20)
   reconstructed_strands = pool.map_async(reconstruct, clusters).get()
   pool.close()
   return reconstructed_strands
def DGSM_Anchored(cluster):
    global f
    f = DoubleSidedGreedyMedian
    return AnchorReconstruct(cluster)
def NW_Anchored(cluster):
    global f
    global ALG
    f = ANCHORED_reconstruct
    ALG = 2
    return AnchorReconstruct(cluster)
def DBMA_Anchored(cluster):
    global f
    global ALG
    f = ANCHORED_reconstruct
    ALG = 0
    return AnchorReconstruct(cluster)
def BMA_Anchored(cluster):
    global f
    global ALG
    f = ANCHORED_reconstruct
    ALG = 1
    return AnchorReconstruct(cluster)
def ORIGINAL_DBMA(cluster):
    global f
    global ALG
    strand_length = STRAND_LENGTHS[0]
    ALG = 0
    return ORIGINAL_reconstruct(cluster)

def g(x):
    if x == 0:
        return DBMA_Anchored
    elif x == 1:
        return BMA_Anchored
    elif x == 2:
        return NW_Anchored
    elif x == 3:
        return DSGM_Anchored

# print(strand_length)
# Recon Algo: use [DBMA BMA NW DGSM]_Anchored
# RECON_ALGO = ORIGINAL_DBMA

parser = argparse.ArgumentParser('DNA consensus')
parser.add_argument('--i', type=str, default="ClusteredStrands.txt", help="input file")
parser.add_argument('--o', type=str, default="ReconstructedStrands.txt", help="output file")
parser.add_argument('--L', type=int, default=120)
parser.add_argument('--W', type=int, default=5)
parser.add_argument('--ALG', type=int, default=0)
parser.add_argument('--E', type=int, default=0)
parser.add_argument('--path', type=str)
parser.add_argument('--coverage', type=int, default=10000)
parser.add_argument('--anchor', type=str)
args = parser.parse_args()
f = None
ALG = None
RECON_ALGO = g(args.ALG)
anchor = args.anchor
# Anchor length
k = len(anchor)
# Search length
D = 14
# l: total strand length including anchor, STRAND_LENGTHS: [w/o anchor, w anchor]
l = args.L
STRAND_LENGTHS = [l, (l-k)//2]
# Strand length of 1 side of the anchor
# strand_length = (strand_length - k)//2
input_file = open(args.i,'r')
data = input_file.readlines()
#clusters = [list(map(lambda x: x.strip(), (data[i*(n+1)+1 : i*(n+1)+S+1]))) for i in range(len(data)//(n+1))]
clusters = []
curr_cluster = []
max_coverage = args.coverage
count = 0
for l in data[2:]:
    # if ((l[1] == '0') OR (l[1] == '1') OR (l[1] == '2') OR (l[1] == '3') OR (l[1] == '4') OR (l[1] == '5') OR (l[1] == '6') OR (l[1] == '7') OR (l[1] == '8') OR (l[1] == '9')):
    if (l[0] in ['0','1','2','3','4','5','6','7','8','9']):
        clusters.append(curr_cluster)
        curr_cluster = []
        count = 0
    elif count >= max_coverage:
        continue
    else:
        curr_cluster.append(l.strip())
        count += 1

clusters.append(curr_cluster)
print("num clusters:", len(clusters))
print("max coverage:", max_coverage)
pool = Pool()
start = time.time()
reconstructed_strands = pool.map(RECON_ALGO, clusters)
print(len(reconstructed_strands[0]))
print("Time:\n\t", time.time() -start)

# Writing Recon strands to file
output_file = open(args.o,'w')
for line in reconstructed_strands:
    output_file.write(line)
    output_file.write('\n')