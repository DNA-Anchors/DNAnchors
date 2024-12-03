
L = ['168']
anchors = ['GTCTGAAC']
for i in range(len(L)):
    file = open('../data/output/baba_smol_wAnchors/EncodedStrands.txt','r')
    all_strands = file.readlines()
    all_strands = list(map(lambda x: x.strip(), all_strands))

    anchor = anchors[i]
    # print('baba_smol_L'+L[i]+'golden added '+anchor)
    output_file = open('../data/output/baba_smol_wAnchors/EncodedStrands_wAnchors.txt', 'w')
    for s in all_strands:
        output_file.write(s[:len(s)//2]+anchor+s[len(s)//2:]+'\n')
