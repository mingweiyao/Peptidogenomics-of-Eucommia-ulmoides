import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def run_cmd(cmd):
    subprocess.run(cmd, check=True)

def parse_outfmt6(path):
    with open(path) as f:
        for line in f:
            yield line.split()

def fasta_id_set(path):
    return {r.id for r in SeqIO.parse(path, "fasta")}

def load_fasta_dict(path, need_ids=False):
    seqs, ids = {}, []
    for r in SeqIO.parse(path, "fasta"):
        seqs[r.id] = str(r.seq)
        if need_ids:
            ids.append(r.id)
    return (seqs, ids) if need_ids else seqs

def make_blast_db(fa, dbtype, out):
    run_cmd(["makeblastdb", "-in", fa, "-dbtype", dbtype, "-out", out, "-parse_seqids"])

def build_hmm(query_fa, hmm="model.hmm", sto="aligned_query.sto"):
    run_cmd(["muscle", "-in", query_fa, "-out", sto, "-moltype", "protein"])
    run_cmd(["hmmbuild", hmm, sto])

def hmmsearch(hmm, protein_fa, out="hmm.tbl", evalue="1e-5"):
    run_cmd(["hmmsearch", "-E", evalue, "--tblout", out, hmm, protein_fa])
    hits = []
    with open(out) as f:
        for l in f:
            if not l.startswith("#"):
                cols = l.split()
                hits.append({ "id": cols[0], 'sequence_id': cols[0], 'start': int(cols[7]), 'end': int(cols[8]), 'strand': cols[9] })
    return hits

def blastp(query_fa, db, out="blastp.tsv", evalue="1e-5"):
    run_cmd(["blastp", "-query", query_fa, "-db", db, "-outfmt", "6", "-evalue", evalue, "-out", out])
    return [{"id": c[1], "sequence_id": c[1], "start": 0, "end": 0, "strand": 0} for c in parse_outfmt6(out)]

def tblastn(query_fa, db, out="tblastn.tsv", evalue="1e-5"):
    run_cmd(["tblastn", "-query", query_fa, "-db", db, "-outfmt", "6", "-evalue", evalue, "-out", out])
    hits = []
    for c in parse_outfmt6(out):
        sseq = c[1]
        s, e = int(c[8]), int(c[9])
        strand = "+" if s < e else "-"
        start, end = (s, e) if s < e else (e, s)
        hits.append({
            "id": f"{sseq}_{start}_{end}_{strand}",
            "sequence_id": sseq,
            "start": start,
            "end": end,
            "strand": strand
        })
    return hits

def remove_duplicates(evidence_sets, chr_id, max_len_aa = 150):
    unique_evidence = []
    in_chr_evidence = [evidence for evidence in evidence_sets if evidence['sequence_id'] in chr_id]
    not_in_chr_evidence = [evidence for evidence in evidence_sets if evidence['sequence_id'] not in chr_id]
    
    unique_not_in_chr_evidence= []
    id_not_chr = set()
    for evidence in not_in_chr_evidence:
        if evidence['id'] not in id_not_chr:
            unique_not_in_chr_evidence.append(evidence)
            id_not_chr.add(evidence['id'])

    groups = defaultdict(list)
    for e in in_chr_evidence:
        groups[(e["sequence_id"], e["strand"])].append(e)
    kept_all = []
    for (seqid, strand), arr in groups.items():
        removed = set()
        for i in range(len(arr)):
            for j in range(i+1, len(arr)):
                a, b = arr[i], arr[j]
                a_s, a_e = sorted((a["start"], a["end"]))
                b_s, b_e = sorted((b["start"], b["end"]))
                if a_e < b_s or b_e < a_s:
                    continue
                if (a_s <= b_s and a_e >= b_e) or (b_s <= a_s and b_e >= a_e):
                    a_len = (a_e - a_s + 1) / 3
                    b_len = (b_e - b_s + 1) / 3
                    if a_len < max_len_aa and b_len < max_len_aa:
                        drop = b if a_len >= b_len else a
                    elif a_len < max_len_aa and b_len >= max_len_aa:
                        drop = b
                    elif b_len < max_len_aa and a_len >= max_len_aa:
                        drop = a
                    else:
                        drop = a if a_len >= b_len else b
                    removed.add(drop["id"])
        kept_all.extend([x for x in arr if x["id"] not in removed])             

    unique_evidence.extend(kept_all)
    unique_evidence.extend(unique_not_in_chr_evidence)
    return unique_evidence

def evidence_to_protein(e, chr_ids, genome, proteins, drop_stop=False):
    sid = e["sequence_id"]
    if sid in chr_ids:
        dna = Seq(genome[sid][e["start"]-1:e["end"]])
        if e["strand"] == "-":
            dna = dna.reverse_complement()
        prot = dna.translate(to_stop=False)
        if drop_stop and "*" in str(prot):
            return None
        return prot
    return Seq(proteins[sid])

def write_hits_fasta(old_fa, evidences, chr_ids, genome, proteins, out="total_hits.fasta", drop_stop=False):
    seqs = {r.id: r.seq for r in SeqIO.parse(old_fa, "fasta")}
    for e in evidences:
        if e["id"] not in seqs:
            p = evidence_to_protein(e, chr_ids, genome, proteins, drop_stop)
            if p:
                seqs[e["id"]] = p
    SeqIO.write([SeqRecord(s, id=i, description="") for i, s in seqs.items()], out, "fasta")

def run_iteration(query_fa, hmm, protein_fa, chr_ids, genome, proteins, max_iter=10, drop_stop=False):
    evidences = []
    prev_ids = fasta_id_set(query_fa)  
    for i in range(max_iter):
        print(f"Starting iteration {i + 1}")
        cur_fa = query_fa if i == 0 else "total_hits.fasta"    
        if i == 0:
            evidences.extend(hmmsearch(hmm, protein_fa))
        evidences.extend(blastp(cur_fa, "protein_db"))
        evidences.extend(tblastn(cur_fa, "genome_db"))
        evidences = remove_duplicates(evidences, chr_ids)
        write_hits_fasta(cur_fa, evidences, chr_ids, genome, proteins, drop_stop=drop_stop)    
        curr_ids = fasta_id_set("total_hits.fasta")
        if curr_ids == prev_ids:
            print("Stop: no new IDs.")
            break
        prev_ids = curr_ids

def main():
    query_fa   = ""
    genome_fa  = ""
    protein_fa = ""
    # 1. 初始化查询序列库
    genome, chr_ids = load_fasta_dict(genome_fa, need_ids=True)
    proteins = load_fasta_dict(protein_fa)
    # 2. 标准化数据库
    make_blast_db(protein_fa, "prot", "protein_db")
    make_blast_db(genome_fa,  "nucl", "genome_db")
    # 3. 构建hmm模型
    build_hmm(query_fa)
    # 4. 开始迭代
    run_iteration(query_fa, "model.hmm", protein_fa, chr_ids, genome, proteins, drop_stop=False)

if __name__ == "__main__":
    main()