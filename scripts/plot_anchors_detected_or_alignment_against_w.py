import math
from tqdm import tqdm
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
    ties = 0
    for i in range(row):
        if matrix[i][-1] > max_score:
            max_score = matrix[i][-1]
            end = i
            ties = 1
        elif matrix[i][-1] == max_score:
            ties += 1
    
    return (ties, max_score)

anchor = 'GTTGACCA'
k = len(anchor)
# Search length
# D = 14
# l: total strand length including anchor, STRAND_LENGTHS: [w/o anchor, w anchor]
l = 176
STRAND_LENGTHS = [l, (l-k)//2]
# Strand length of 1 side of the anchor
# strand_length = (strand_length - k)//2
input_file = open('../data/output/baba_smol_wAnchors_P0.02_N10/UnderlyingClusters.txt','r')
data = input_file.readlines()
#clusters = [list(map(lambda x: x.strip(), (data[i*(n+1)+1 : i*(n+1)+S+1]))) for i in range(len(data)//(n+1))]
clusters = []
curr_cluster = []
max_coverage = math.inf
count = 0
for l in data[1:]:
    # if ((l[1] == '0') OR (l[1] == '1') OR (l[1] == '2') OR (l[1] == '3') OR (l[1] == '4') OR (l[1] == '5') OR (l[1] == '6') OR (l[1] == '7') OR (l[1] == '8') OR (l[1] == '9')):
    if (l[1]=='L'):
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

import matplotlib.pyplot as plt
ties = {}
max_scores = {}
x = []
y = []
for D in tqdm(range(8,176,2)):
    ties = {}
    max_scores = {}
    for c in clusters:
        for s in c:
            n = len(s)
            if D <= n:
                middle = s[(n//2)-(D//2):(n//2)-(D//2)+D]
            else:
                middle = s
            t,m = NW(middle, anchor)
            if t not in ties:
                ties[t] = 1
            else:
                ties[t]+=1
            if m not in max_scores:
                max_scores[m]=1
            else:
                max_scores[m]+=1
    x.append(D)
    # Number of anchors detected
    if 1 in ties:
        y.append(ties[1]/2550)
    else:
        y.append(0)

    # # Avg Score of Anchor
    # t = 0
    # for p in max_scores:
    #     t+= p*max_scores[p]
    # y.append(t/2550)

plt.xlabel("Window Size", fontsize=16)
plt.ylabel("Distinct Anchors Detected (%)", fontsize=16)
# plt.ylabel('Avg Alignment Score', fontsize=16)
plt.plot(x,y)
plt.savefig("Distinct_vs_w.pdf", bbox_inches='tight')
# plt.savefig("Alignment_Score_vs_w.pdf", bbox_inches='tight')
print('peak at: ', x[y.index(max(y))])
