from itertools import product

ALL_FILES = ['../data/output/baba_smol_wAnchors/EncodedStrands.txt']

k = 8 #anchor length
D = 168

all_kmers = set()
for name in ALL_FILES:
    file = open(name,'r')
    all_strands = file.readlines()
    print("num strands:",len(all_strands))
    strand_len = len(all_strands[0].strip())
    start = (strand_len//2)- (D//2)
    end = start+ D
    all_strands = list(map(lambda x: x.strip()[start:end], all_strands))
    for s in all_strands:
        for i in range(D-k+1):
            all_kmers.add(s[i:i+k])
            

    

nucleotides = ['A', 'T', 'C', 'G']
entire_space = {''.join(p) for p in product(nucleotides, repeat=k)}

# print('k:',k)
# print('D:',D)


ALL_CANDIDATE_ANCHORS = set()
for s in (entire_space-all_kmers):
    a,c,g,t = s.count('A'),s.count('C'),s.count('G'),s.count('T')
    if a == c == g == t:
        ALL_CANDIDATE_ANCHORS.add(s)

print("Anchors with balanced characters:\n\t", len(ALL_CANDIDATE_ANCHORS))

# print(ALL_CANDIDATE_ANCHORS)