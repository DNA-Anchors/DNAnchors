import configparser
import sys
from reedsolo import RSCodec, ReedSolomonError
import time
import random
from multiprocessing import Pool

BASE_DIR=sys.argv[1]
config_location = sys.argv[2]
config = configparser.ConfigParser()
config.read(config_location)
encode_or_decode  = int(sys.argv[3])  # 0: encode, 1: decode
out_dir=sys.argv[4]
skip_ECC  = int(sys.argv[5])  # 0: encode with RS, 1: encode without RS
parameters = config['parameters']
num_rows = int(parameters['num_rows'])
symbol_size  = int(parameters['symbol_size'])
code_length =  255
fec_l=int(parameters['redundancy'])
data_length = code_length-fec_l
file_size = int(parameters['file_size']) #if specified for encoding, uses only the first {file_size} bytes form the file
#mapping=int(parameters['mapping_scheme'])

file_locations = config['file_locations']
inputfile=(file_locations['input_f'])
outputfile_encoded=(file_locations['output_f_encoded'])
reconstructed_strands=(file_locations['reconstructed_f'])
outputfile_decoded=(file_locations['output_f_decoded'])
priority_file=(file_locations['priority_f'])


print("Symbol size: ",str(symbol_size)+" bits")
print("Using: ","("+str(code_length)+","+str(fec_l)+") ReedSolomon Codes")
# Num bases used for indexing
ind = 4
strand_length = ind + ((symbol_size//2)*num_rows)
print("Payload Length: "+str(strand_length))


# Encode
start = time.time()
if not encode_or_decode:
    RS_matrix = [[] for i in range(num_rows)]
    rsc = RSCodec(fec_l)
    with open(BASE_DIR+'/input/'+inputfile, 'rb') as file:
        curr_row = 0
        while True:
            data = bytearray(file.read(data_length*(symbol_size//8)))
            if not data:
                break
            elif len(data) < data_length*(symbol_size//8):
                data += b'\x00' * (data_length*(symbol_size//8)-len(data))
            random.seed(curr_row)
            XOR_mask = bytearray([random.randint(0,255) for _ in range(len(data))])
            RS_matrix[curr_row] = rsc.encode(bytearray(a ^ b for a, b in zip(data, XOR_mask)))
            curr_row += 1
    for r in range(curr_row, num_rows):
        random.seed(r)
        data = b'\x00' * data_length*(symbol_size//8)
        XOR_mask = bytearray([random.randint(0,255) for _ in range(len(data))])
        RS_matrix[r] = rsc.encode(bytearray(a ^ b for a, b in zip(data, XOR_mask)))
    file.close()
    ENCODING_LOOKUP_TABLE = {0: 'AAAA', 1: 'AAAC', 2: 'AAAG', 3: 'AAAT', 4: 'AACA', 5: 'AACC', 6: 'AACG', 7: 'AACT', 8: 'AAGA', 9: 'AAGC', 10: 'AAGG', 11: 'AAGT', 12: 'AATA', 13: 'AATC', 14: 'AATG', 15: 'AATT', 16: 'ACAA', 17: 'ACAC', 18: 'ACAG', 19: 'ACAT', 20: 'ACCA', 21: 'ACCC', 22: 'ACCG', 23: 'ACCT', 24: 'ACGA', 25: 'ACGC', 26: 'ACGG', 27: 'ACGT', 28: 'ACTA', 29: 'ACTC', 30: 'ACTG', 31: 'ACTT', 32: 'AGAA', 33: 'AGAC', 34: 'AGAG', 35: 'AGAT', 36: 'AGCA', 37: 'AGCC', 38: 'AGCG', 39: 'AGCT', 40: 'AGGA', 41: 'AGGC', 42: 'AGGG', 43: 'AGGT', 44: 'AGTA', 45: 'AGTC', 46: 'AGTG', 47: 'AGTT', 48: 'ATAA', 49: 'ATAC', 50: 'ATAG', 51: 'ATAT', 52: 'ATCA', 53: 'ATCC', 54: 'ATCG', 55: 'ATCT', 56: 'ATGA', 57: 'ATGC', 58: 'ATGG', 59: 'ATGT', 60: 'ATTA', 61: 'ATTC', 62: 'ATTG', 63: 'ATTT', 64: 'CAAA', 65: 'CAAC', 66: 'CAAG', 67: 'CAAT', 68: 'CACA', 69: 'CACC', 70: 'CACG', 71: 'CACT', 72: 'CAGA', 73: 'CAGC', 74: 'CAGG', 75: 'CAGT', 76: 'CATA', 77: 'CATC', 78: 'CATG', 79: 'CATT', 80: 'CCAA', 81: 'CCAC', 82: 'CCAG', 83: 'CCAT', 84: 'CCCA', 85: 'CCCC', 86: 'CCCG', 87: 'CCCT', 88: 'CCGA', 89: 'CCGC', 90: 'CCGG', 91: 'CCGT', 92: 'CCTA', 93: 'CCTC', 94: 'CCTG', 95: 'CCTT', 96: 'CGAA', 97: 'CGAC', 98: 'CGAG', 99: 'CGAT', 100: 'CGCA', 101: 'CGCC', 102: 'CGCG', 103: 'CGCT', 104: 'CGGA', 105: 'CGGC', 106: 'CGGG', 107: 'CGGT', 108: 'CGTA', 109: 'CGTC', 110: 'CGTG', 111: 'CGTT', 112: 'CTAA', 113: 'CTAC', 114: 'CTAG', 115: 'CTAT', 116: 'CTCA', 117: 'CTCC', 118: 'CTCG', 119: 'CTCT', 120: 'CTGA', 121: 'CTGC', 122: 'CTGG', 123: 'CTGT', 124: 'CTTA', 125: 'CTTC', 126: 'CTTG', 127: 'CTTT', 128: 'GAAA', 129: 'GAAC', 130: 'GAAG', 131: 'GAAT', 132: 'GACA', 133: 'GACC', 134: 'GACG', 135: 'GACT', 136: 'GAGA', 137: 'GAGC', 138: 'GAGG', 139: 'GAGT', 140: 'GATA', 141: 'GATC', 142: 'GATG', 143: 'GATT', 144: 'GCAA', 145: 'GCAC', 146: 'GCAG', 147: 'GCAT', 148: 'GCCA', 149: 'GCCC', 150: 'GCCG', 151: 'GCCT', 152: 'GCGA', 153: 'GCGC', 154: 'GCGG', 155: 'GCGT', 156: 'GCTA', 157: 'GCTC', 158: 'GCTG', 159: 'GCTT', 160: 'GGAA', 161: 'GGAC', 162: 'GGAG', 163: 'GGAT', 164: 'GGCA', 165: 'GGCC', 166: 'GGCG', 167: 'GGCT', 168: 'GGGA', 169: 'GGGC', 170: 'GGGG', 171: 'GGGT', 172: 'GGTA', 173: 'GGTC', 174: 'GGTG', 175: 'GGTT', 176: 'GTAA', 177: 'GTAC', 178: 'GTAG', 179: 'GTAT', 180: 'GTCA', 181: 'GTCC', 182: 'GTCG', 183: 'GTCT', 184: 'GTGA', 185: 'GTGC', 186: 'GTGG', 187: 'GTGT', 188: 'GTTA', 189: 'GTTC', 190: 'GTTG', 191: 'GTTT', 192: 'TAAA', 193: 'TAAC', 194: 'TAAG', 195: 'TAAT', 196: 'TACA', 197: 'TACC', 198: 'TACG', 199: 'TACT', 200: 'TAGA', 201: 'TAGC', 202: 'TAGG', 203: 'TAGT', 204: 'TATA', 205: 'TATC', 206: 'TATG', 207: 'TATT', 208: 'TCAA', 209: 'TCAC', 210: 'TCAG', 211: 'TCAT', 212: 'TCCA', 213: 'TCCC', 214: 'TCCG', 215: 'TCCT', 216: 'TCGA', 217: 'TCGC', 218: 'TCGG', 219: 'TCGT', 220: 'TCTA', 221: 'TCTC', 222: 'TCTG', 223: 'TCTT', 224: 'TGAA', 225: 'TGAC', 226: 'TGAG', 227: 'TGAT', 228: 'TGCA', 229: 'TGCC', 230: 'TGCG', 231: 'TGCT', 232: 'TGGA', 233: 'TGGC', 234: 'TGGG', 235: 'TGGT', 236: 'TGTA', 237: 'TGTC', 238: 'TGTG', 239: 'TGTT', 240: 'TTAA', 241: 'TTAC', 242: 'TTAG', 243: 'TTAT', 244: 'TTCA', 245: 'TTCC', 246: 'TTCG', 247: 'TTCT', 248: 'TTGA', 249: 'TTGC', 250: 'TTGG', 251: 'TTGT', 252: 'TTTA', 253: 'TTTC', 254: 'TTTG', 255: 'TTTT'}
    def get_DNA_strand_8bit(col):
        strand = [ENCODING_LOOKUP_TABLE[r[col]] for r in RS_matrix]
        return ''.join(strand)
    def get_DNA_strand_16bit(col):
        strand = [ENCODING_LOOKUP_TABLE[r[2*col]]+ENCODING_LOOKUP_TABLE[r[2*col+1]] for r in RS_matrix]
        return ''.join(strand)
    pool = Pool()
    if symbol_size == 8:
        encoded_strands = pool.map(get_DNA_strand_8bit, range(code_length))
    else:
        encoded_strands = pool.map(get_DNA_strand_16bit, range(code_length))

    with open(BASE_DIR+'/output/'+out_dir+'/'+outputfile_encoded, 'w') as file:
        for i in range(len(encoded_strands)):
            file.write(ENCODING_LOOKUP_TABLE[i] + encoded_strands[i] + '\n')
    file.close()

# Decode
else:
    rsc = RSCodec(fec_l)
    DECODING_LOOKUP_TABLE = {'AAAA': 0, 'AAAC': 1, 'AAAG': 2, 'AAAT': 3, 'AACA': 4, 'AACC': 5, 'AACG': 6, 'AACT': 7, 'AAGA': 8, 'AAGC': 9, 'AAGG': 10, 'AAGT': 11, 'AATA': 12, 'AATC': 13, 'AATG': 14, 'AATT': 15, 'ACAA': 16, 'ACAC': 17, 'ACAG': 18, 'ACAT': 19, 'ACCA': 20, 'ACCC': 21, 'ACCG': 22, 'ACCT': 23, 'ACGA': 24, 'ACGC': 25, 'ACGG': 26, 'ACGT': 27, 'ACTA': 28, 'ACTC': 29, 'ACTG': 30, 'ACTT': 31, 'AGAA': 32, 'AGAC': 33, 'AGAG': 34, 'AGAT': 35, 'AGCA': 36, 'AGCC': 37, 'AGCG': 38, 'AGCT': 39, 'AGGA': 40, 'AGGC': 41, 'AGGG': 42, 'AGGT': 43, 'AGTA': 44, 'AGTC': 45, 'AGTG': 46, 'AGTT': 47, 'ATAA': 48, 'ATAC': 49, 'ATAG': 50, 'ATAT': 51, 'ATCA': 52, 'ATCC': 53, 'ATCG': 54, 'ATCT': 55, 'ATGA': 56, 'ATGC': 57, 'ATGG': 58, 'ATGT': 59, 'ATTA': 60, 'ATTC': 61, 'ATTG': 62, 'ATTT': 63, 'CAAA': 64, 'CAAC': 65, 'CAAG': 66, 'CAAT': 67, 'CACA': 68, 'CACC': 69, 'CACG': 70, 'CACT': 71, 'CAGA': 72, 'CAGC': 73, 'CAGG': 74, 'CAGT': 75, 'CATA': 76, 'CATC': 77, 'CATG': 78, 'CATT': 79, 'CCAA': 80, 'CCAC': 81, 'CCAG': 82, 'CCAT': 83, 'CCCA': 84, 'CCCC': 85, 'CCCG': 86, 'CCCT': 87, 'CCGA': 88, 'CCGC': 89, 'CCGG': 90, 'CCGT': 91, 'CCTA': 92, 'CCTC': 93, 'CCTG': 94, 'CCTT': 95, 'CGAA': 96, 'CGAC': 97, 'CGAG': 98, 'CGAT': 99, 'CGCA': 100, 'CGCC': 101, 'CGCG': 102, 'CGCT': 103, 'CGGA': 104, 'CGGC': 105, 'CGGG': 106, 'CGGT': 107, 'CGTA': 108, 'CGTC': 109, 'CGTG': 110, 'CGTT': 111, 'CTAA': 112, 'CTAC': 113, 'CTAG': 114, 'CTAT': 115, 'CTCA': 116, 'CTCC': 117, 'CTCG': 118, 'CTCT': 119, 'CTGA': 120, 'CTGC': 121, 'CTGG': 122, 'CTGT': 123, 'CTTA': 124, 'CTTC': 125, 'CTTG': 126, 'CTTT': 127, 'GAAA': 128, 'GAAC': 129, 'GAAG': 130, 'GAAT': 131, 'GACA': 132, 'GACC': 133, 'GACG': 134, 'GACT': 135, 'GAGA': 136, 'GAGC': 137, 'GAGG': 138, 'GAGT': 139, 'GATA': 140, 'GATC': 141, 'GATG': 142, 'GATT': 143, 'GCAA': 144, 'GCAC': 145, 'GCAG': 146, 'GCAT': 147, 'GCCA': 148, 'GCCC': 149, 'GCCG': 150, 'GCCT': 151, 'GCGA': 152, 'GCGC': 153, 'GCGG': 154, 'GCGT': 155, 'GCTA': 156, 'GCTC': 157, 'GCTG': 158, 'GCTT': 159, 'GGAA': 160, 'GGAC': 161, 'GGAG': 162, 'GGAT': 163, 'GGCA': 164, 'GGCC': 165, 'GGCG': 166, 'GGCT': 167, 'GGGA': 168, 'GGGC': 169, 'GGGG': 170, 'GGGT': 171, 'GGTA': 172, 'GGTC': 173, 'GGTG': 174, 'GGTT': 175, 'GTAA': 176, 'GTAC': 177, 'GTAG': 178, 'GTAT': 179, 'GTCA': 180, 'GTCC': 181, 'GTCG': 182, 'GTCT': 183, 'GTGA': 184, 'GTGC': 185, 'GTGG': 186, 'GTGT': 187, 'GTTA': 188, 'GTTC': 189, 'GTTG': 190, 'GTTT': 191, 'TAAA': 192, 'TAAC': 193, 'TAAG': 194, 'TAAT': 195, 'TACA': 196, 'TACC': 197, 'TACG': 198, 'TACT': 199, 'TAGA': 200, 'TAGC': 201, 'TAGG': 202, 'TAGT': 203, 'TATA': 204, 'TATC': 205, 'TATG': 206, 'TATT': 207, 'TCAA': 208, 'TCAC': 209, 'TCAG': 210, 'TCAT': 211, 'TCCA': 212, 'TCCC': 213, 'TCCG': 214, 'TCCT': 215, 'TCGA': 216, 'TCGC': 217, 'TCGG': 218, 'TCGT': 219, 'TCTA': 220, 'TCTC': 221, 'TCTG': 222, 'TCTT': 223, 'TGAA': 224, 'TGAC': 225, 'TGAG': 226, 'TGAT': 227, 'TGCA': 228, 'TGCC': 229, 'TGCG': 230, 'TGCT': 231, 'TGGA': 232, 'TGGC': 233, 'TGGG': 234, 'TGGT': 235, 'TGTA': 236, 'TGTC': 237, 'TGTG': 238, 'TGTT': 239, 'TTAA': 240, 'TTAC': 241, 'TTAG': 242, 'TTAT': 243, 'TTCA': 244, 'TTCC': 245, 'TTCG': 246, 'TTCT': 247, 'TTGA': 248, 'TTGC': 249, 'TTGG': 250, 'TTGT': 251, 'TTTA': 252, 'TTTC': 253, 'TTTG': 254, 'TTTT': 255}
    RS_matrix_T = [None for i in range(code_length)]
    with open(BASE_DIR+'/output/'+out_dir+'/'+reconstructed_strands,'r') as file:
        # all_strands = map(lambda x: x.strip(), file.readlines())
        
        #below for simulating missing m strands
        m = 1 
        all_strands = list(map(lambda x: x.strip(), file.readlines()))[m:]
        for s in all_strands:
            RS_matrix_T[DECODING_LOOKUP_TABLE[s[:ind]]] = s[ind:]
    file.close()
    successfully_decoded = [0]*num_rows
    num_erasures = 0
    for i in range(len(RS_matrix_T)):
        if RS_matrix_T[i] == None:
            num_erasures += 1
            RS_matrix_T[i] = 'AAAA'*num_rows*(symbol_size//8)
    def get_symbols_8bit(col):
        symbols = []
        for i in range(code_length):
            symbols.append(DECODING_LOOKUP_TABLE[RS_matrix_T[i][col*4:col*4+4]])
        symbols = bytearray(symbols)
        try:
            result = rsc.decode(symbols)
            data = result[0]
            successfully_decoded[col] = len(result[2])
        except:
            data = symbols
            successfully_decoded[col] = 'False'
        random.seed(col)
        XOR_mask = bytearray([random.randint(0,255) for _ in range(len(data))])
        return bytearray(a ^ b for a, b in zip(data, XOR_mask))
    def get_symbols_16bit(col):
        symbols = []
        for i in range(code_length):
            symbols.append(DECODING_LOOKUP_TABLE[RS_matrix_T[i][col*8:col*8+4]])
            symbols.append(DECODING_LOOKUP_TABLE[RS_matrix_T[i][col*8+4:col*8+8]])
        symbols = bytearray(symbols)
        try: 
            result = rsc.decode(bytearray(symbols))
            data = result[0]
            successfully_decoded[col] = len(result[2])
        except:
            data = symbols
            successfully_decoded[col] = 'False'
        random.seed(col)
        XOR_mask = bytearray([random.randint(0,255) for _ in range(len(data))])
        return bytearray(a ^ b for a, b in zip(data, XOR_mask))
    # pool = Pool()
    # print('CODE LENGTH',code_length)
    decoded_rows = []
    if symbol_size == 8:
        # decoded_rows = pool.map(get_symbols_8bit, range(num_rows))
        for i in range(num_rows):
            decoded_rows.append(get_symbols_8bit(i))
    else:
        decoded_rows = pool.map(get_symbols_16bit, range(num_rows))


    num_complete_rows = file_size//(data_length*(symbol_size//8))
    partial_row = file_size%(data_length*(symbol_size//8))
    with open(BASE_DIR+'/output/'+out_dir+'/'+outputfile_decoded, 'wb') as file:
        for i in range(num_complete_rows):
            file.write(decoded_rows[i])
        if partial_row:
            file.write(decoded_rows[num_complete_rows][:partial_row])
        # for r in decoded_rows:
        #     file.write(r)
    file.close()

    with open(BASE_DIR+'/output/'+out_dir+'/'+'decoder_log.txt','w') as file:
        file.write("Number of Erasures:\t" + str(num_erasures)+'\n')
        data = [[i,'True' if type(successfully_decoded[i])!=str else successfully_decoded[i],successfully_decoded[i] if type(successfully_decoded[i])==int else '-'] for i in range(num_rows)]
        data = [[" Row ", " Decoded? ", " Num Errors "]] + data

        col_widths = [max(len(str(item)) for item in column) for column in zip(*data)]
        file.write('+' + '+'.join('-' * (col_widths[i] + 2) for i in range(len(data[0]))) + '+\n')
        file.write('|' + ' | '.join(f"{item:<{col_widths[i]}}" for i, item in enumerate(data[0])) + '|\n')
        file.write('+' + '+'.join('-' * (col_widths[i] + 2) for i in range(len(data[0]))) + '+\n')
        for row in data[1:]:
            file.write('|' + ' | '.join(f"{item:<{col_widths[i]}}" for i, item in enumerate(row)) + '|\n')
        file.write('+' + '+'.join('-' * (col_widths[i] + 2) for i in range(len(data[0]))) + '+\n')
    file.close()

print("Time:\n\t", time.time() -start)
