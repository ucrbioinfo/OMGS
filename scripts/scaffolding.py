import sys

def read_keyfile(myfile2, sname_qry):
    key_file = myfile2 + "_key.txt"
    with open(key_file) as f_key:
        for i in range(0, 4): # 4 header lines
            f_key.readline()  
        for line in f_key:
            line = line.strip()
            cols = line.split('\t')
            qry_id = int(cols[0])
            seq_name = cols[1]
            sname_qry[seq_name] = qry_id
    f_key.close()

def revcomp(seq):
    rcseq = ""
    for i in range(0,len(seq)):
        c = seq[len(seq)-1-i]
        if c=='A':
            rcseq = rcseq + 'T'
        elif c=='T':
            rcseq = rcseq + 'A'
        elif c=='G':
            rcseq = rcseq + 'C'
        elif c=='C':
            rcseq = rcseq + 'G'
        elif c=='a':
            rcseq = rcseq + 't'
        elif c=='t':
            rcseq = rcseq + 'a'
        elif c=='g':
            rcseq = rcseq + 'c'
        elif c=='c':
            rcseq = rcseq + 'g'
        else:
            rcseq = rcseq + c 
    return rcseq


def scaffolding(myfile2_list, fasta_file, orders_file, gaps_file, ori_file, output_dir):
    print '---------------cutting seqs-------------------'
    # read key_file
#    sname_qry = {}
#    for i in range(len(myfile2_list)):
#        myfile2 = myfile2_list[i]
#        read_keyfile(myfile2, sname_qry)

    # read fasta file
    seqs_list = {}
    num = 1
    for line in open(fasta_file):
        line = line.strip()
        if line[0] == '>':
            seq_name = line[1:]
#            current_name = seq_name
            current_id = num
            num += 1
            seqs_list[current_id] = []
        else:
            seqs_list[current_id].append(line)
    
    seqs = {}
    for qry in seqs_list:
        seq_list = seqs_list[qry]
        whole_seq = ''.join(seq_list)
        seqs[qry] = whole_seq

    # read orders
    used_qry = set([])
    orders = []
    for line in open(orders_file):
        line = line.strip()
        cols = line.split()
        order = []

        for item in cols:
            used_qry.add(int(item))
            order.append(int(item))
        orders.append(order)

    # read gaps
    gaps = {}
    for line in open(gaps_file):
        line = line.strip()
        cols = line.split()
        qry1 = int(cols[0])
        qry2 = int(cols[1])
        gap = float(cols[2])
        gaps[qry1, qry2] = gap

    # read oris
    oris = {}
    for line in open(ori_file):
        line = line.strip()
        cols = line.split()
        qry = int(cols[0])
        ori = cols[1]
        oris[qry] = ori
        

    scaffolds = []
    for order in orders:
        scaffold = ""
        for i in range(len(order)):
            qry = order[i]
            seq = seqs[qry]
            if oris[qry] == '-':
                seq = revcomp(seq)
            if scaffold == "":
                scaffold += seq
            else:
                last_qry = order[i-1]
                gap = int(gaps[last_qry, qry])
                if gap >= 0:
                    scaffold += 'N'*gap 
                else:
                    scaffold += 'N'*100
#                scaffold += "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
                scaffold += seq
        scaffolds.append(scaffold)
    for qry in seqs:
        if qry not in used_qry:
            seq = seqs[qry]
            scaffolds.append(seq)

    fo = file(output_dir+"/scaffolds.fasta", 'w')
    n = 1
    for seq in scaffolds:
        fo.write(">scaffold_"+str(n)+"\n")
        n+=1
        fo.write(seq+"\n")



