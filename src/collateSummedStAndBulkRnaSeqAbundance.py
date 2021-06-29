import re

mylog = open(snakemake.log.log,'wt')
sys.stdout=mylog
sys.stderr=mylog

db = {}
st = open(snakemake.input.st, 'rt')
head = {}
valset = set()
sumst = 0
for row in st:
    row = row.strip("\n").split("\t")
    if len(db) == 0:
        for n,val in enumerate(row):
            if val != "":
                id = val
                if snakemake.params.doGeneLevel:
                    id = re.sub("_i\d+\Z", "", id)
                if id not in db.keys():
                    db[id] = { "st" : 0 }
                    head[id] = []
                elif val in valset:
                    print("Error: transcript {} appears twice in st RNAseq results".format(val))
                    sys.exit(-1)
                valset.add(val)
                head[id].append(n)
    else:
        for id in db.keys():
            for val in head[id]:
                db[id]["st"] += float(row[val])
                sumst += float(row[val])
st.close()

bulk = open(snakemake.input.bulk, 'rt')
valset = set()
sumbulk = 0
k=0
for row in bulk:
    row = row.strip().split("\t")
    val = row[0]
    id = val
    if snakemake.params.doGeneLevel:
        id = re.sub("_i\d+\Z", "", id)
    if id not in db.keys():
        print("Warning: id {i}, with bulk abundance {a}, is missing in st results".format(i=id, a=row[1]))
        db[id] = { "st" : -1, "bulk" : 0 }
        k += 1
    elif "bulk" not in db[id].keys():
        db[id]["bulk"] = 0
    db[id]["bulk"] += float(row[1])
    sumbulk += float(row[1])
    if val in valset:
        print("Error: transcript {} appears twice in bulk RNAseq results". format(id))
        sys.exit(-1)
    valset.add(val)
bulk.close()

out = open(snakemake.output.tsv, 'wt')
out.write("{}\tcount_st\tcount_bulk\tCPM_st\tCPM_bulk\n".format("genes" if snakemake.params.doGeneLevel else "transcripts"))

n = 0
m = 0
for id in db.keys():
    if "bulk" not in db[id].keys():
        print("Warning: transcript {i}, with st abundance {a}, is missing in bulk RNAseq results".format(i=id, a = db[id]["st"]))
        db[id]["bulk"] = -1
        m+=1
    elif "st" in db[id].keys():
        n+=1
    out.write("{id}\t{cst}\t{cb}\t{tst}\t{tb}\n".format(id=id,
                                                        cst=db[id]["st"],
                                                        cb=db[id]["bulk"],
                                                        tst="{:.2E}".format(db[id]["st"]/sumst * 1E6if db[id]["st"] != -1 else -1),
                                                        tb="{:.2E}".format(db[id]["bulk"]/sumbulk * 1E6 if db[id]["bulk"] != -1 else -1)))
print("Found {n} common {feature} between st and bulk seq".format(n=n, feature="genes" if snakemake.params.doGeneLevel else "transcripts"))
print("{n} {feature} was missing from st RNAseq".format(n=k, feature="genes" if snakemake.params.doGeneLevel else "transcripts"))
print("{n} {feature} was missing from bulk RNAseq".format(n=m, feature="genes" if snakemake.params.doGeneLevel else "transcripts"))

out.close()
mylog.close()

