import hashlib
import base64
from struct import *


# Sha1coin block 432100 - BlockHash : 32da761289f5602eed03971a0ae398eded5a4a6b69d5c3be0c2bff48de2da16a
hdr_version = 2
hdr_hashPrevBlock = "c21b38bf0c08c269a3a01e1f4509ae3101d1ee32beb06272739f2b48c071b335"
hdr_hashMerkleRoot = "40ea7f4acd9bfeca05ca48244cd76712662ec239820c5dfec5fedebe16f2ca50"
hdr_time = 1456070206
hdr_bits = "1c521ddf"
hdr_nonce = 736338799


################################################################################

data = pack('L', hdr_version)
data += hdr_hashPrevBlock.decode('hex')[::-1]
data += hdr_hashMerkleRoot.decode('hex')[::-1]
data += pack('L', hdr_time)
data += hdr_bits.decode('hex')[::-1]
data += pack('L', hdr_nonce)


# SHA-256 double (for blockhash)

hash = hashlib.sha256(hashlib.sha256(data).digest()).digest()
print "BlockHash:    " + hash[::-1].encode('hex')


# sha1coinhash (for mining)

str = base64.b64encode(hashlib.sha1(data).digest())[:26]
str += str

hash256 = 0

for i in range(26):
    subhash = hashlib.sha1(str[i:i+12]).digest()
    subhash256 = 0
    
    for j in range(20):
        subhash256 += int(subhash[j].encode('hex'), 16) << (j+12)*8
        
    hash256 ^= subhash256
    
hash_str = ""

for i in range(32):
    hash_str += '%02x' % (hash256 >> ((31-i)*8) & 0xff)

print "MiningHash:   " + hash_str


# Target

target256 = int(hdr_bits[2:], 16) * 2 ** (8 * (int(hdr_bits[:2], 16) - 3))

target256_str = ""

for i in range(32):
    target256_str += '%02x' % (target256 >> ((31-i)*8) & 0xff)

print "Target:       " + target256_str


