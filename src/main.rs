use num_bigint::{BigInt, BigUint};
use num_traits::Num;
use sha2::{Digest, Sha256};
use crate::elliptic::{EllipticCryptosystem, EllipticCurve};

mod elliptic;
mod primes;

fn main() {
    let mut crypto = EllipticCryptosystem::new_standard("P-256").unwrap();
    let private_key = crypto.generate_private_key();
    let message = "Hello, world!";
    let digest: [u8; 32] = Sha256::digest(message.as_bytes()).into();
    let signature = crypto.sign(&digest, &private_key);
    let verification_result = crypto.verify(&digest, &signature, &private_key.get_public_key()).unwrap();
    println!("verification result: {}", verification_result);
    
    let other_private_key = crypto.generate_private_key();
    let other_verification_result = crypto.verify(&digest, &signature, &other_private_key.get_public_key()).unwrap();
    println!("verification result: {}", other_verification_result);
    
}
