import struct, sys
def rint(f): return struct.unpack("<i", f.read(4))[0]
def rstr(f):
    n=0; s=0
    while True:
        b=f.read(1)[0]; n|=(b&0x7f)<<s; s+=7
        if not b&0x80: break
    return f.read(n).decode("latin-1")
f=open(sys.argv[1],"rb")
assert f.read(3)==b"BPM"; f.read(1); ver=rint(f)
if ver & 0x1000: ver^=0x1000
rstr(f)
if ver>1: rstr(f)
n=rint(f); f.read(4*n)
names=[rstr(f) for _ in range(n)]
f.read(n)
out=open(sys.argv[2],"w"); out.write("name\tchrom\tpos\tsnp\n")
for i in range(n):
    v=rint(f); rstr(f); name=rstr(f)
    for _ in range(3): rstr(f)
    f.read(4)
    for _ in range(2): rstr(f)
    snp=rstr(f); chrom=rstr(f)
    for _ in range(2): rstr(f)
    pos=rstr(f)
    for _ in range(2): rstr(f)
    f.read(8)
    for _ in range(7): rstr(f)
    f.read(3); f.read(1)
    if v>=7: f.read(16)
    if v>=8: rstr(f)
    out.write(f"{name}\t{chrom}\t{pos}\t{snp}\n")
print("loci",n,"version",ver,"last locus version",v)
