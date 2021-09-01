from tensorflow import keras
import numpy as np
import csv
import argparse

def one_hot(seq):
    one_hot_seq = np.zeros((4, len(seq)))
    for i in range(len(seq)):
            if seq[i] == "A" or seq[i] == "a":
                one_hot_seq[0][i] = 1
            if seq[i] == "C" or seq[i] == "c":
                one_hot_seq[1][i] = 1
            if seq[i] == "G" or seq[i] == "g":
                one_hot_seq[2][i] = 1
            if seq[i] == "T" or seq[i] == "t":
                one_hot_seq[3][i] = 1
    return one_hot_seq

def make_prediction(seq, i, k):
    one_hot_seq = one_hot(seq[i:i+k])
    mid = one_hot(seq[i+int(k/2)])
    one_hot_seq = np.expand_dims(one_hot_seq, axis=-1)
    one_hot_seq = np.expand_dims(one_hot_seq, axis=0)
    mf_pred = np.exp(model_mf.predict(one_hot_seq)[0][0])
    sub_pred = model_sub.predict(one_hot_seq)
    sub_pred = np.multiply(1-mid[:,0],sub_pred[0])
    sub_pred = sub_pred/np.sum(sub_pred)
    return mf_pred, sub_pred


parser = argparse.ArgumentParser(description='Predicts mutation frequency and substitution profile.')
parser.add_argument('input', metavar='input', type=str,
                    help='path to the input fasta file')
parser.add_argument('-o', metavar='output', type=str, default='output.csv',
                    help='path to the output csv file')
parser.add_argument('-k', metavar='k', type=int, default=15, choices=[5, 9, 15, 21], 
                    help='k-mer size')

args = parser.parse_args()

header = ['sequence id', 'k-mer start', 'k-mer end', 'middle nucleotide position', 
            'middle nucleotide', 'mutation frequency', 'N>A', 'N>C', 'N>G', 'N>T']
with open(args.o, 'w') as f: 
    # using csv.writer method from CSV package
    write = csv.writer(f)
    write.writerow(header)

k = args.k
model_mf = keras.models.load_model('./models/model_'+str(k)+'_mf.h5')
model_sub = keras.models.load_model('./models/model_'+str(k)+'_sub.h5')

new_seq = 0
seq = ''
inputs = 0
results=[]
seq_id = 0
outputs = []
with open(args.input,'r') as data:
    for line in csv.reader(data):
        if line[0][0] == '>':
            if seq != '':
                if len(seq) >= k:
                    for i in range(len(seq) - k+1):
                        mf_pred, sub_pred = make_prediction(seq, i, k)
                        kmer_data = [seq_id, i+1, i+k, i+int(k/2)+1, seq[i+int(k/2)], mf_pred, sub_pred[0], sub_pred[1], sub_pred[2], sub_pred[3]]
                        with open(args.o, 'a') as f: 
                            write = csv.writer(f)
                            write.writerow(kmer_data)
            new_seq = 1
            seq_id += 1
        elif line[0][0] in ['A', 'C', 'G', 'T']:
            if new_seq == 1:
                seq = line[0]
                new_seq = 0
            else:
                seq += line[0]

if len(seq) >= k:
    for i in range(len(seq) - k+1):
        mf_pred, sub_pred = make_prediction(seq, i, k)
        kmer_data = [seq_id, i+1, i+k, i+int(k/2)+1, seq[i+int(k/2)], mf_pred, sub_pred[0], sub_pred[1], sub_pred[2], sub_pred[3]]
        with open(args.o, 'a') as f: 
            write = csv.writer(f)
            write.writerow(kmer_data)
