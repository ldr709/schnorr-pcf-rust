use arrayvec::ArrayVec;
use core::mem::swap;
use curve25519_dalek::constants::ED25519_BASEPOINT_TABLE;
use curve25519_dalek::edwards::{CompressedEdwardsY, EdwardsPoint};
use curve25519_dalek::scalar::Scalar;
use ed25519_dalek::{SignatureError, VerifyingKey};
use rand::{CryptoRng, RngCore};
use sha2::{Sha512, digest::Digest};

use crate::dl_pcf::*;

#[derive(Copy, Clone)]
pub struct PubKey(CompressedEdwardsY);

pub struct Signer {
    prover: Prover,
    verifier: Verifier,
    pk: PubKey,
    sk_i: Scalar,
}

pub struct SignerMsg(Proof);

#[allow(non_snake_case)]
pub struct SignerState<'a> {
    signer: &'a Signer,
    msg: &'a [u8],
    r_i: Scalar,
    R_i: EdwardsPoint,
}

#[allow(non_snake_case)]
#[derive(Copy, Clone)]
pub struct Signature {
    s: Scalar,
    R: CompressedEdwardsY,
}

pub struct SignatureShare(Signature);

impl PubKey {
    #[allow(non_snake_case)]
    pub fn verify(&self, msg: &[u8], sig: &Signature) -> bool {
        let mut hasher = Sha512::new();
        hasher.update(sig.R.as_bytes());
        hasher.update(self.0.as_bytes());
        hasher.update(msg);
        let h = Scalar::from_bytes_mod_order_wide(&hasher.finalize().into());
        let R_verify = EdwardsPoint::vartime_double_scalar_mul_basepoint(
            &-h, &self.0.decompress().unwrap(), &sig.s);

        R_verify.compress() == sig.R
    }
}

impl TryFrom<PubKey> for VerifyingKey {
    type Error = SignatureError;
    fn try_from(pk: PubKey) -> Result<Self, Self::Error> {
        VerifyingKey::from_bytes(&pk.0.to_bytes())
    }
}

impl From<Signature> for ed25519_dalek::Signature {
    fn from(sig: Signature) -> Self {
        ed25519_dalek::Signature::from_components(sig.R.to_bytes(), sig.s.to_bytes())
    }
}

impl Signer {
    #[allow(non_snake_case)]
    pub fn round1<'a>(&'a self, msg: &'a [u8]) -> (SignerMsg, SignerState<'a>) {
        let (r_i, R_i, proof) = self.prover.prove(msg);
        (SignerMsg(proof), SignerState { signer: self, msg, r_i, R_i })
    }
}

impl<'a> SignerState<'a> {
    #[allow(non_snake_case)]
    pub fn round2(&self, signer_msg: SignerMsg) -> Result<SignatureShare, InvalidProof> {
        let R_j = self.signer.verifier.verify(self.msg, signer_msg.0)?;
        let R = (self.R_i + R_j).compress();

        let mut hasher = Sha512::new();
        hasher.update(R.as_bytes());
        hasher.update(self.signer.pk.0.as_bytes());
        hasher.update(self.msg);
        let h = Scalar::from_bytes_mod_order_wide(&hasher.finalize().into());

        let s_i = h * self.signer.sk_i + self.r_i;
        Ok(SignatureShare(Signature { s: s_i, R }))
    }
}

pub fn keygen<R: RngCore + CryptoRng>(rng: &mut R) -> (PubKey, [Signer; 2]) {
    let sk = [Scalar::random(rng), Scalar::random(rng)];
    let pk = &(sk[0] + sk[1]) * ED25519_BASEPOINT_TABLE;
    let pk = PubKey(pk.compress());

    let ell = 3072;
    let mut prover_verifier_pairs = [setup(rng, ell), setup(rng, ell)];

    // Swap verifier halves, to give each party half of each prover-verifier pair.
    let pvps_split = prover_verifier_pairs.split_at_mut(1);
    swap(&mut pvps_split.0[0].1, &mut pvps_split.1[0].1);

    let arr_vec: ArrayVec<_, 2> = prover_verifier_pairs.into_iter().enumerate().map(|(i, pv)| {
        Signer {
            prover: pv.0,
            verifier: pv.1,
            pk,
            sk_i: sk[i],
        }
    }).collect();
    (pk, arr_vec.into_inner().ok().unwrap())
}

pub fn combine_signature_shares(shares: &[SignatureShare; 2]) -> Signature {
    assert_eq!(shares[0].0.R, shares[1].0.R);
    Signature {
        s: shares[0].0.s + shares[1].0.s,
        R: shares[0].0.R,
    }
}

pub fn sign_2party(parties: [Signer; 2], msg: &[u8]) -> Signature {
    let (round1_msgs, round1_states): (ArrayVec<_, 2>, ArrayVec<_, 2>) = parties.iter().map(|p| {
        p.round1(msg)
    }).unzip();
    let sig_shares: ArrayVec<_, 2> = round1_states.iter().zip(round1_msgs.into_iter().rev())
    .map(|(s, m)| {
        s.round2(m).unwrap()
    }).collect();
    combine_signature_shares(&sig_shares.into_inner().ok().unwrap())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_chacha::ChaCha20Rng;
    use rand::SeedableRng;

    #[test]
    fn signing() {
        let mut rng = ChaCha20Rng::from_entropy();
        let (pk, signers) = keygen(&mut rng);
        println!("Sampled key.");
        let message = b"Test message";
        let signature = sign_2party(signers, message);
        println!("Computed signature.");

        assert!(pk.verify(message, &signature));

        let dalek_pk = VerifyingKey::try_from(pk).unwrap();
        dalek_pk.verify_strict(message, &signature.into()).unwrap();
    }
}
