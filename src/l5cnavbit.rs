//! # GPS L5 CNAV navigation messages (thin wrapper over `CNavBit`)
//!
//! GPS L5 CNAV carries the *same* CNAV message payload (MT10 / MT11 / MT30 …) as L2C —
//! IS-GPS-705 and IS-GPS-200 differ only in the transmission channel (PRN, symbol rate),
//! not in the message content or FEC. `CNavBit::get_frame_data` already produces both:
//! `param != 0` selects the L5 timing (6 s messages) and the L5 convolutional-encoder
//! continuity state.
//!
//! The previous standalone implementation re-derived the message bit packing and the FEC,
//! and diverged from the standard in several places — wrong Keplerian scale factors and
//! field widths (audit H9 / H10 / H11), `sqrtA` encoded instead of ΔA, and a wrong
//! convolutional polynomial. It is now a thin wrapper that delegates to the verified
//! `CNavBit` encoder; `nav_decode::decode_cnav_ephemeris` round-trips that encoder.

use crate::cnavbit::CNavBit;
use crate::types::GnssTime;
use crate::{GpsAlmanac, GpsEphemeris, IonoParam, UtcParam};

/// GPS L5 CNAV navigation-bit generator. Delegates to [`CNavBit`] with the L5 `param`.
#[derive(Clone)]
pub struct L5CNavBit {
    inner: CNavBit,
}

impl Default for L5CNavBit {
    fn default() -> Self {
        Self::new()
    }
}

impl L5CNavBit {
    pub fn new() -> Self {
        Self {
            inner: CNavBit::new(),
        }
    }

    /// Generates one L5 CNAV message (600 FEC symbols). `param` is ignored — L5 always uses
    /// the `param = 1` path of `CNavBit` (6 s message cadence + L5 convolutional state).
    pub fn get_frame_data(
        &mut self,
        start_time: GnssTime,
        svid: i32,
        _param: i32,
        nav_bits: &mut [i32; 600],
    ) -> i32 {
        self.inner.get_frame_data(start_time, svid, 1, nav_bits)
    }

    pub fn set_ephemeris(&mut self, svid: i32, eph: &GpsEphemeris) -> i32 {
        self.inner.set_ephemeris(svid, eph)
    }

    pub fn set_almanac(&mut self, alm: &[GpsAlmanac; 32]) -> i32 {
        self.inner.set_almanac(alm)
    }

    pub fn set_iono_utc(&mut self, iono_param: &IonoParam, utc_param: &UtcParam) -> i32 {
        self.inner.set_iono_utc(iono_param, utc_param)
    }
}
