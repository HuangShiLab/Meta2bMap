// Shared 64-bit tag/k-mer hash type.
pub type Hash = u64;

// FNV-1a hash for canonical byte-string tags.
pub fn hash_bytes(bytes: &[u8]) -> Hash {
    let mut hash: u64 = 0xcbf29ce484222325; // FNV offset basis
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x100000001b3); // FNV prime
    }
    hash
}
