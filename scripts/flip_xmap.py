def flip_xmap(prefix):
    in_file = prefix + ".xmap"
    out_file = prefix + "_flip.xmap"
    fout = open(out_file, 'w')
    for line in open(in_file):
        if line[0] == '#':
            fout.write(line)
            continue
        line = line.strip()
        cols = line.split("\t")
        sign = cols[7]
        if sign == '+':
            q_start = cols[5]
            q_end = cols[6]
            r_start = cols[3]
            r_end = cols[4]
        else:
            q_start = cols[6]
            q_end = cols[5]
            r_start = cols[4]
            r_end = cols[3]
        outline = cols[0] + "\t" + cols[2] + "\t" + cols[1] + "\t" + q_start + "\t" + q_end + "\t" + r_start + "\t" + r_end + "\t" + cols[7] + "\t" + cols[8] + "\t" + cols[9] + "\t" + cols[11] + "\t" + cols[10] + "\t"+ cols[12] + "\t" + cols[13] + "\n"
        fout.write(outline)



