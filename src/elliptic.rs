use std::cmp::min;
use crate::primes;
use num::Integer;
use num_bigint::{BigInt, BigUint, RandBigInt};
use num_traits::{Num, One, Signed, Zero};
use sha2::{Digest, Sha256};
use std::ops::Add;
use std::rc::Rc;
use rand::CryptoRng;

#[derive(thiserror::Error, Debug)]
pub enum EllipticError {
    #[error("the signature is malformed")]
    MalformedSignature,
    #[error("point is not on the curve")]
    NotOnCurve,
    #[error("unexpected error: {0}")]
    #[allow(dead_code)]
    Other(String),
}

pub trait EllipticCurve: PartialEq + Sized + Clone {
    fn has_point(&self, p: &Point) -> bool;
    fn add_points(&self, p1: &CurvePoint<Self>, p2: &CurvePoint<Self>) -> CurvePoint<Self>;
    fn multiply_point(&self, n: &BigInt, point: &CurvePoint<Self>) -> CurvePoint<Self> {
        if n.is_zero() {
            panic!("cannot multiply point by zero without knowing the generator")
        }
        // adjust for negative multiple
        let (mut n, mut point): (BigUint, CurvePoint<Self>) = if n.is_negative() {
            (BigUint::try_from(n.abs()).unwrap(), self.negate(point))
        } else {
            (BigUint::try_from(n).unwrap(), point.clone())
        };

        let mut result = point.clone();

        let zero = BigUint::zero();
        let one = BigUint::one();
        let two = &one + &one;
        // since we start with result = point (or -point) instead of neutral element of a group
        // (which we cannot represent anyway as it's at infinity), we subtract 1 from exponent
        n -= one;

        while n != zero {
            if n.is_odd() {
                result = &result + &point;
            }
            n = &n / &two;
            point = &point + &point;
        }
        result
    }
    fn negate(&self, p1: &CurvePoint<Self>) -> CurvePoint<Self>;
}

struct WeierstrassParametersInternal {
    p: BigUint,
    a: BigUint,
    b: BigUint,
    // for fast comparison between curves in different forms
    curve_hash: [u8; 32],
}
#[derive(Clone)]
pub struct WeierstrassCurve {
    parameters: Rc<WeierstrassParametersInternal>,
}

#[derive(Debug, Clone)]
pub struct Point {
    pub x: BigUint,
    pub y: BigUint,
}

#[derive(Debug, Clone)]
pub struct CurvePoint<C: EllipticCurve> {
    pub x: BigUint,
    pub y: BigUint,
    curve: C,
}

impl<C: EllipticCurve> PartialEq for &CurvePoint<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.curve != other.curve {
            panic!("comparing points on different curves")
        }
        self.x == other.x && self.y == other.y
    }
}

impl<C: EllipticCurve> Add for &CurvePoint<C> {
    type Output = CurvePoint<C>;
    fn add(self, other: &CurvePoint<C>) -> CurvePoint<C> {
        self.curve.add_points(&self, &other)
    }
}
impl PartialEq for WeierstrassCurve {
    fn eq(&self, other: &Self) -> bool {
        self.parameters.curve_hash == other.parameters.curve_hash
    }
}

impl EllipticCurve for WeierstrassCurve {
    fn has_point(&self, point: &Point) -> bool {
        let p = &self.parameters.p;
        let a = &self.parameters.a;
        let b = &self.parameters.b;
        let x = &point.x;
        let y = &point.y;

        // no point in using modpow here due to small exponent
        let lhs = (y * y) % p;
        let rhs = ((x * x * x) % p + (a * x) % p + b % p) % p;
        lhs == rhs
    }

    fn add_points(&self, p1: &CurvePoint<Self>, p2: &CurvePoint<Self>) -> CurvePoint<Self> {
        self.add_points(p1, p2)
    }

    fn negate(&self, p: &CurvePoint<Self>) -> CurvePoint<Self> {
        CurvePoint {
            x: p.x.clone(),
            y: &self.parameters.p - p.y.clone(),
            curve: self.clone(),
        }
    }
}

impl WeierstrassCurve {
    pub fn new(p: BigUint, a: BigUint, b: BigUint) -> WeierstrassCurve {
        if !primes::Verification::is_prime(&p) {
            panic!("p is not prime")
        }

        let mut hasher = Sha256::new();
        hasher.update(b"WeierstrassCurve");
        hasher.update(Sha256::digest(a.to_bytes_be()));
        hasher.update(Sha256::digest(b.to_bytes_be()));
        let result = hasher.finalize();
        let parameters = Rc::new(WeierstrassParametersInternal {
            p,
            a,
            b,
            curve_hash: result
                .as_slice()
                .try_into()
                .expect("unexpected hash length"),
        });
        WeierstrassCurve { parameters }
    }

    fn from_string_parameters(p_str: &str, a_str: &str, b_str: &str) -> Self {
        let p = BigUint::from_str_radix(p_str, 16).unwrap();
        // just a sanity check
        if !primes::Verification::is_prime(&p) {
            panic!("p is not prime")
        }
        let a = BigUint::from_str_radix(a_str, 16).unwrap();
        let b = BigUint::from_str_radix(b_str, 16).unwrap();
        WeierstrassCurve::new(p, a, b)
    }
    pub fn new_standard_curve(name: &str) -> Result<Self, EllipticError> {
        match name {
            "P-256" => Ok(WeierstrassCurve::from_string_parameters(
                "ffffffff00000001000000000000000000000000ffffffffffffffffffffffff",
                "ffffffff00000001000000000000000000000000fffffffffffffffffffffffc",
                "5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b",
            )),
            _ => Err(EllipticError::Other(format!("unknown curve: {}", name))),
        }
    }
    fn has_curve_point(&self, p: &CurvePoint<Self>) -> bool {
        // we ensure that the point is on curve when we create it, so only thing we need to check
        // is whether we are actually using the same curve
        &p.curve == self
    }

    pub fn get_point(&self, p: Point) -> Result<CurvePoint<Self>, EllipticError> {
        if !self.has_point(&p) {
            return Err(EllipticError::NotOnCurve);
        }
        Ok(CurvePoint {
            x: p.x,
            y: p.y,
            curve: self.clone(),
        })
    }

    pub fn double_point(&self, point: &CurvePoint<Self>) -> CurvePoint<Self> {
        let p = &self.parameters.p;
        let a = &self.parameters.a;
        let b = &self.parameters.b;

        let x = &point.x;
        let y = &point.y;

        let one = BigUint::one();
        let two = &one + &one;
        let three = &two + &one;

        let slope = (((three * x * x % p) + a) % p) * (&two * y).modinv(p).unwrap();
        let new_x = (&slope * &slope + (p - &two * x % p)) % p;
        let new_y = (&slope * (x + p - &new_x) + p - y) % p;
        CurvePoint {
            x: new_x,
            y: new_y,
            curve: self.clone(),
        }
    }
    fn add_points(&self, p1: &CurvePoint<Self>, p2: &CurvePoint<Self>) -> CurvePoint<Self> {
        if p1 == p2 {
            return self.double_point(p1);
        }
        let p = &self.parameters.p;
        let a = &self.parameters.a;
        let b = &self.parameters.b;

        let x1 = &p1.x;
        let y1 = &p1.y;
        let x2 = &p2.x;
        let y2 = &p2.y;
        let slope = ((y2 + p - y1) * (x2 + p - x1).modinv(p).unwrap()) % p;
        let new_x = (&slope * &slope + (p - x1) + (p - x2)) % p;
        let new_y = (&slope * (x1 + p - &new_x) - y1) % p;
        CurvePoint {
            x: new_x,
            y: new_y,
            curve: self.clone(),
        }
    }
}

#[derive(Clone)]
pub struct EllipticCryptosystem<C: EllipticCurve> {
    rng: rand::rngs::ThreadRng, // ThreadRng is crypto sage
    curve: C,
    generator: CurvePoint<C>,
    order: BigUint,
}

impl EllipticCryptosystem<WeierstrassCurve> {
    pub fn new_standard(name: &str) -> Result<Self, EllipticError> {
        let curve = WeierstrassCurve::new_standard_curve(name)?;
        let (generator, order) = match name {
            "P-256" => (
                curve.get_point(Point {
                    x: BigUint::from_str_radix(
                        "6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296",
                        16,
                    )
                    .unwrap(),
                    y: BigUint::from_str_radix(
                        "4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5",
                        16,
                    )
                    .unwrap(),
                }).unwrap(),
                BigUint::from_str_radix(
                    "ffffffff00000000ffffffffffffffffbce6faada7179e84f3b9cac2fc632551",
                    16,
                )
                .unwrap(),
            ),
            _ => {
                return Err(EllipticError::Other(format!(
                    "unknown standard cryptosystem: {}",
                    name
                )))
            }
        };
        Ok(EllipticCryptosystem {
            rng: rand::thread_rng(),
            curve,
            generator,
            order,
        })
    }
}

#[derive(Clone)]
pub struct ECPrivateKey<C: EllipticCurve> {
    private_exponent: BigUint,
    public_key: ECPublicKey<C>,
}

#[derive(Clone)]
pub struct ECPublicKey<C: EllipticCurve> {
    public_point: CurvePoint<C>,
}


impl<C: EllipticCurve> ECPrivateKey<C> {
    pub fn new(cryptosystem: &EllipticCryptosystem<C>, private_exponent: BigUint) -> Self {
        let public_point = cryptosystem.curve.multiply_point(&BigInt::from(private_exponent.clone()), &cryptosystem.generator);
        ECPrivateKey { private_exponent, public_key: ECPublicKey { public_point} }
    }
    
    pub fn get_public_key(&self) -> ECPublicKey<C> {
        self.public_key.clone()
    }
}

#[derive(Clone)]
pub struct ECDSASignature {
    r: BigUint,
    s: BigUint,
}
impl<C: EllipticCurve> EllipticCryptosystem<C> {
    pub fn generate_private_key(&mut self) -> ECPrivateKey<C> {
        let d = self.rng.gen_biguint_below(&self.order);
        ECPrivateKey::new(self, d)
    }

    //pub fn load_public_key()
    pub fn sign(&mut self, hash: &[u8], private_key: &ECPrivateKey<C>) -> ECDSASignature {
        // L_n is the bit length of the group order n
        let ln = self.order.bits() as usize;
        if ln % 8 != 0 {
            panic!("unsupported generator order in the cryptosystem; must be divisible by 8");
        }
        let ln = ln / 8;
        let ln_capped = std::cmp::min(ln, hash.len());
        // Z are the L_n leftmost bits of the hash.
        let z = BigUint::from_bytes_be(&hash[..ln_capped]);
        loop {
            let k = self.rng.gen_biguint_below(&self.order);
            let signature_point = self.curve.multiply_point(&BigInt::from(k.clone()), &self.generator);
            let CurvePoint { x: x1, y: y1, .. } = signature_point;
            let r = x1 % &self.order;
            if r.is_zero() { continue; } // this is astronomically unlikely
            let s = k.modinv(&self.order).unwrap()*(&z + &r*&private_key.private_exponent) % &self.order;
            return ECDSASignature { r, s };
        }
    }

    pub fn verify(&mut self, hash: &[u8], signature: &ECDSASignature, public_key: &ECPublicKey<C>) -> Result<bool, EllipticError> {
        let r = &signature.r;
        let s = &signature.s;
        if r.is_zero() || r >= &self.order || s.is_zero() || s >= &self.order {
            return Err(EllipticError::MalformedSignature);
        }
        let s_inv = match s.modinv(&self.order) {
            Some(inv) => inv,
            None => return Err(EllipticError::MalformedSignature),
        };
        // L_n is the bit length of the group order n
        let ln = self.order.bits() as usize;
        if ln % 8 != 0 {
            panic!("unsupported generator order in the cryptosystem; must be divisible by 8");
        }
        let ln = ln / 8;
        let ln_capped = std::cmp::min(ln, hash.len());
        // Z are the L_n leftmost bits of the hash.
        let z = BigUint::from_bytes_be(&hash[..ln_capped]);

        let u1 = (z * &s_inv) % &self.order;
        let u2 = (r * &s_inv) % &self.order;
        let verification_point = (
            &self.curve.multiply_point(&BigInt::from(u1), &self.generator) +
            &self.curve.multiply_point(&BigInt::from(u2), &public_key.public_point)
        );
        Ok(r == &verification_point.x)
    }
}