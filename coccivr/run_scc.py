# SPDX-License-Identifier: GPL-3.0-or-later
import os
import subprocess

from .configure import exe


def ensure_blast_db(db):
    """Build BLAST indices next to the FASTA if they are not there yet."""
    if os.path.exists(db + '.nin'):
        return
    proc = subprocess.run([exe['makeblastdb'], '-dbtype', 'nucl', '-in', db, '-out', db],
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            'Could not build BLAST database for {0} (is the directory writable?):\n{1}'
            .format(db, proc.stderr.strip()))


def run_scc_blastn(query,db,cpus,n_hits=3000,topOnly=True,ident=90, cov=60,sum_cov=60):
    ensure_blast_db(db)
    cmd = [exe['blastn'],
           '-num_threads', str(cpus),
           '-task', 'blastn',
           '-query', query,
           '-db', db,
           '-max_target_seqs', str(n_hits),
           '-outfmt', '6 qaccver saccver pident length mismatch gapopen '
                      'qstart qend sstart send evalue score qlen slen']
    bsn = []
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout :
        part = line.strip().split('\t')
        if len(part) < 14 :
            continue
        part[2:] = [float(v) for v in part[2:]]
        if part[2] < ident :
            continue
        part[5] = '+'
        if part[9] < part[8] :
            part[8:10], part[5] = part[9:7:-1], '-'
        elif part[7] < part[6] :
            part[6:8], part[5] = part[7:5:-1], '-'

        if part[9]-part[8]+1 < 30:
            continue
        if (part[9]-part[8]+1)*100. < cov*part[13] :
            continue
        bsn.append(part)
    p.communicate()

    if topOnly :
        bsn = sorted(bsn, key=lambda x:[x[0], x[6], x[7]]) 
        
        for i1, b1 in enumerate(bsn[:-1]) :
            for b2 in bsn[i1+1:] :
                if b2[0] == '' :
                    continue
                if b1[0] != b2[0] or b1[7] < b2[6] :
                    break
                o = min(b1[7], b2[7]) - max(b1[6], b2[6]) + 1
                if o >= 0.4*(b1[7]-b1[6]+1) or o > 0.4*(b2[7]-b2[6]+1) :
                    if b1[2] >= b2[2]+1 and b1[11] >= b2[11] :
                        b2[0] = ''
                    elif b2[2] >= b1[2]+1 and b2[11] >= b1[11] :
                        b1[0] = ''
        bsn = sorted([b for  b in bsn if b[0] != ''], key=lambda x:(x[1], x[8], x[9]))
        

        scores = {}
        for b in bsn :
            if b[1] not in scores :
                scores[b[1]] = [b[-1], [b[8], b[9]]]
            elif b[8] - 1 > scores[b[1]][-1][1] :
                scores[b[1]].append([b[8], b[9]])
            elif scores[b[1]][-1][1] < b[9] :
                scores[b[1]][-1][1] = b[9]
        for r, s in scores.items() :
            s[1:] = [sum([ss[1]-ss[0]+1 for ss in s[1:]])]
        bsn = [ b for b in bsn if scores[b[1]][1]*100 >= scores[b[1]][0]*sum_cov ]

        bsn = sorted(bsn, key=lambda x: (x[0], x[6], x[7]))
        for i1, b1 in enumerate(bsn[:-1]):
            for b2 in bsn[i1 + 1:]:
                if b2[0] == '':
                    continue
                if b1[0] != b2[0] or b1[7] < b2[6]:
                    break
                o = min(b1[7], b2[7]) - max(b1[6], b2[6]) + 1
                if o >= 0.4 * (b1[7] - b1[6] + 1) or o > 0.4 * (b2[7] - b2[6] + 1):
                    if (b1[2], b1[11]) > (b2[2], b2[11]):
                        b2[0] = ''
                    else :
                        b1[0] = ''


        bsn = sorted([b for b in bsn if b[0] != ''], key=lambda x: (x[1], x[8], x[9]))
        outputs = {}
        for b in bsn:
            if b[1] not in outputs:
                outputs[b[1]] = [b[-1], [b[2]*(b[9]-b[8]+1), (b[9]-b[8]+1)], [[b[8], b[9]]], [b]]
            else :
                outputs[b[1]][3].append(b)
                outputs[b[1]][1][0] += b[2]*(b[9]-b[8]+1)
                outputs[b[1]][1][1] += (b[9]-b[8]+1)
                if b[8] - 1 > outputs[b[1]][2][-1][1]:
                    outputs[b[1]][2].append([b[8], b[9]])
                elif outputs[b[1]][2][-1][1] < b[9]:
                    outputs[b[1]][2][-1][1] = b[9]
        for r in list(outputs.keys()) :
            match = outputs[r]
            match[1] = match[1][0]/match[1][1]
            match[2] = sum([ m[1]-m[0]+1 for m in match[2] ])
            if match[2]*100 < match[0]*sum_cov :
                outputs.pop(r)
        #print(outputs)
    
        
        predict_scc_type = []
        for output in outputs:
            scc_tpye = str(output).strip().split('|')[0]
            predict_scc_type.append(scc_tpye)
        scc_type_print = ','.join(predict_scc_type)
        return scc_type_print,outputs

