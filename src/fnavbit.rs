//! # Galileo F/NAV Навигационные сообщения
//!
//! Этот модуль реализует синтез навигационных битов для Galileo F/NAV (Freely accessible Navigation).
//!
//! ## Назначение
//! F/NAV - это формат навигационных сообщений системы Galileo, передаваемый на частоте E5a
//! (1176.45 МГц). Предоставляет свободно доступную навигационную информацию для гражданских
//! пользователей с улучшенными характеристиками по сравнению с GPS.
//!
//! ## Основные функции модуля
//! - Генерация страниц F/NAV с эфемеридной и часовой информацией
//! - Формирование альманаха Galileo для всех спутников созвездия
//! - Кодирование параметров ионосферы и UTC для Galileo System Time (GST)
//! - Формирование сообщений о состоянии сигналов (Signal-in-Space Status)
//! - Поддержка сообщений поиска и спасения (Search and Rescue)
//!
//! F/NAV использует структуру страниц по 244 символа каждая, передаваемых со скоростью
//! 25 символов в секунду, что обеспечивает улучшенную помехоустойчивость и точность.

//----------------------------------------------------------------------
// fnavbit.rs:
//   Implementation of navigation bit synthesis class for F/NAV
//
//          Copyright (C) 2020-2029 by Jun Mo, All rights reserved.
//
//----------------------------------------------------------------------

use crate::types::*;
use crate::COMPOSE_BITS;

// Count number of set bits in a byte
fn count1(n: u8) -> i32 {
    let mut n = n;
    let mut count = 0;
    while n > 0 {
        n &= n - 1;
        count += 1;
    }
    count & 1
}

const SQRT_A0: f64 = 5_440.588_203_494_177;
const NOMINAL_I0: f64 = 0.977_384_381_116_824_6;

#[derive(Clone)]
pub struct FNavBit {
    pub gal_eph_data: [[[u32; 7]; 4]; 36],    // 36 satellites, 4 page types, 7 words each
    pub gal_alm_data: [[[u32; 7]; 2]; 12],    // 12 almanac groups, 2 page types, 7 words each
    pub gal_utc_data: [u32; 4],               // UTC parameters
    pub gal_iono_data: [u32; 2],              // Ionosphere parameters
}

impl FNavBit {
    const SYNC_PATTERN: [i32; 12] = [1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0];

    pub fn new() -> Self {
        FNavBit {
            gal_eph_data: [[[0; 7]; 4]; 36],
            gal_alm_data: [[[0; 7]; 2]; 12],
            gal_utc_data: [0; 4],
            gal_iono_data: [0; 2],
        }
    }

    pub fn get_frame_data(&self, start_time: GnssTime, svid: i32, _param: i32, nav_bits: &mut [i32]) -> i32 {
        // First determine the current TOW and subframe number
        let week = start_time.Week + start_time.MilliSeconds / 604800000;
        let milliseconds = start_time.MilliSeconds % 604800000;
        let tow = milliseconds / 1000;
        let gst = (((week - 1024) & 0xfff) << 20) + tow;
        let subframe = (tow % 1200) / 50; // two round of 600s frame (24 subframes) to hold 36 almanacs
        let page = (tow % 50) / 10;
        
        let mut encode_data = [0u32; 7];
        self.get_page_data(svid, page, subframe, gst as u32, &mut encode_data);
        let crc_result = Self::crc24q_encode(&encode_data, 248);

        // Place message bits and CRC into a single array (248 + 24 = 272 bits)
        let mut uncoded_bits = [0i32; 272];
        for i in 0..248 {
            uncoded_bits[i] = ((encode_data[i / 32] >> (31 - (i % 32))) & 1) as i32;
        }
        for i in 0..24 {
            uncoded_bits[248 + i] = ((crc_result >> (23 - i)) & 1) as i32;
        }

        // FEC (Forward Error Correction) Rate 1/2 Convolutional Encoder for F/NAV
        // G1 = 1110101b (0x75), G2 = 1011011b (0x5B)
        let mut encoded_symbols = [0i32; 544]; // 272 * 2 = 544 after convolutional encoding
        let mut conv_state = 0u8;
        
        for i in 0..272 {
            conv_state = (conv_state << 1) | (uncoded_bits[i] as u8);
            encoded_symbols[i * 2] = count1(conv_state & 0x75);
            encoded_symbols[i * 2 + 1] = count1(conv_state & 0x5B);
        }

        // Add 12-bit sync pattern
        for i in 0..12 {
            nav_bits[i] = Self::SYNC_PATTERN[i];
        }

        // Block interleaving (8 rows, 68 columns)
        for i in 0..544 {
            if 12 + i < nav_bits.len() {
                nav_bits[i + 12] = encoded_symbols[(i % 8) * 68 + (i / 8)];
            }
        }

        0
    }

    pub fn set_ephemeris(&mut self, svid: i32, eph: &GpsEphemeris) -> i32 {
        if !(1..=36).contains(&svid) || eph.valid == 0 {
            return 0;
        }
        let eph_data = &mut self.gal_eph_data[(svid - 1) as usize];
        FNavBit::compose_eph_words(eph, eph_data);
        svid
    }

    pub fn set_almanac(&mut self, alm: &[GpsAlmanac]) -> i32 {
        let mut week = 0i32;
        
        for i in 0..36.min(alm.len()) {
            if (alm[i].valid & 1) != 0 {
                week = alm[i].week;
                break;
            }
        }
        
        for i in 0..12 {
            let start_idx = i * 3;
            if start_idx + 2 < alm.len() {
                FNavBit::compose_alm_words(&alm[start_idx..start_idx + 3], &mut self.gal_alm_data[i], week);
            }
        }
        
        0
    }

    pub fn set_iono_utc(&mut self, iono_param: &IonoParam, utc_param: &UtcParam) -> i32 {
        // Put ionosphere parameters (assuming NeQuick model)
        let uint_value = Self::unscale_uint(iono_param.a0, -2);
        self.gal_iono_data[0] = COMPOSE_BITS!(uint_value, 5, 11);
        
        let int_value = Self::unscale_int(iono_param.a1, -8);
        self.gal_iono_data[0] |= COMPOSE_BITS!(int_value >> 6, 0, 5);
        self.gal_iono_data[1] = COMPOSE_BITS!(int_value, 26, 6);
        
        let int_value = Self::unscale_int(iono_param.a2, -15);
        self.gal_iono_data[1] |= COMPOSE_BITS!(int_value, 12, 14);
        self.gal_iono_data[1] |= COMPOSE_BITS!(iono_param.flag, 7, 5);

        // Put UTC parameters
        let int_value = Self::unscale_int(utc_param.A0, -30);
        self.gal_utc_data[0] = COMPOSE_BITS!(int_value >> 26, 0, 6);
        self.gal_utc_data[1] = COMPOSE_BITS!(int_value, 6, 26);
        
        let int_value = Self::unscale_int(utc_param.A1, -50);
        self.gal_utc_data[1] |= COMPOSE_BITS!(int_value >> 18, 0, 6);
        self.gal_utc_data[2] = COMPOSE_BITS!(int_value, 14, 18);
        self.gal_utc_data[2] |= COMPOSE_BITS!(utc_param.TLS, 6, 8);
        self.gal_utc_data[2] |= COMPOSE_BITS!(utc_param.tot >> 2, 0, 6);
        self.gal_utc_data[3] = COMPOSE_BITS!(utc_param.tot, 30, 2);
        self.gal_utc_data[3] |= COMPOSE_BITS!(utc_param.WN, 22, 8);
        self.gal_utc_data[3] |= COMPOSE_BITS!(utc_param.WNLSF, 14, 8);
        self.gal_utc_data[3] |= COMPOSE_BITS!(utc_param.DN, 11, 3);
        self.gal_utc_data[3] |= COMPOSE_BITS!(utc_param.TLSF, 3, 8);

        0
    }

    fn compose_eph_words(ephemeris: &GpsEphemeris, eph_data: &mut [[u32; 7]; 4]) -> i32 {
        // Initialize all data to zero
        for i in 0..4 {
            for j in 0..7 {
                eph_data[i][j] = 0;
            }
        }

        // PageType 1
        // Note: TOW will be added in get_page_data
        eph_data[0][0] = COMPOSE_BITS!(1, 16, 6); // Page type = 1
        
        let uint_value = ephemeris.toc / 60;
        eph_data[0][1] = COMPOSE_BITS!(ephemeris.iodc, 16, 10) | COMPOSE_BITS!(uint_value, 2, 14);
        
        let int_value = Self::unscale_int(ephemeris.af2, -59);
        eph_data[0][1] |= COMPOSE_BITS!(int_value >> 4, 0, 2);
        eph_data[0][2] = COMPOSE_BITS!(int_value, 28, 4);
        
        let int_value = Self::unscale_int(ephemeris.af1, -46);
        eph_data[0][2] |= COMPOSE_BITS!(int_value, 7, 21);
        
        eph_data[0][3] = COMPOSE_BITS!(ephemeris.health & 0x3, 29, 2); // E5a HS
        let int_value = Self::unscale_int(ephemeris.af0, -34);
        eph_data[0][3] |= COMPOSE_BITS!(int_value >> 2, 0, 29);
        eph_data[0][4] = COMPOSE_BITS!(int_value, 30, 2);
        
        let int_value = Self::unscale_int(ephemeris.tgd, -32);
        eph_data[0][5] = COMPOSE_BITS!(int_value, 11, 10) | COMPOSE_BITS!((ephemeris.health >> 2) & 0x1, 10, 1);
        eph_data[0][6] = 0; // GST will be added in get_page_data

        // PageType 2
        eph_data[1][0] = (2 << 16) | COMPOSE_BITS!(ephemeris.iodc, 6, 10);
        let int_value = Self::unscale_int(ephemeris.M0 / std::f64::consts::PI, -31);
        eph_data[1][0] |= COMPOSE_BITS!(int_value >> 26, 0, 6);
        eph_data[1][1] = COMPOSE_BITS!(int_value, 6, 26);
        
        let int_value = Self::unscale_int(ephemeris.omega_dot / std::f64::consts::PI, -43);
        eph_data[1][1] |= COMPOSE_BITS!(int_value >> 18, 0, 6);
        eph_data[1][2] = COMPOSE_BITS!(int_value, 14, 18);
        
        let uint_value = Self::unscale_uint(ephemeris.ecc, -33);
        eph_data[1][2] |= COMPOSE_BITS!(uint_value >> 18, 0, 14);
        eph_data[1][3] = COMPOSE_BITS!(uint_value, 14, 18);
        
        let uint_value = Self::unscale_uint(ephemeris.sqrtA, -19);
        eph_data[1][3] |= COMPOSE_BITS!(uint_value >> 18, 0, 14);
        eph_data[1][4] = COMPOSE_BITS!(uint_value, 14, 18);
        
        let int_value = Self::unscale_int(ephemeris.omega0 / std::f64::consts::PI, -31);
        eph_data[1][4] |= COMPOSE_BITS!(int_value >> 18, 0, 14);
        eph_data[1][5] = COMPOSE_BITS!(int_value, 14, 18);
        
        let int_value = Self::unscale_int(ephemeris.idot / std::f64::consts::PI, -43);
        eph_data[1][5] |= COMPOSE_BITS!(int_value, 0, 14);
        eph_data[1][6] = 0; // for GST

        // PageType 3
        eph_data[2][0] = (3 << 16) | COMPOSE_BITS!(ephemeris.iodc, 6, 10);
        let int_value = Self::unscale_int(ephemeris.i0 / std::f64::consts::PI, -31);
        eph_data[2][0] |= COMPOSE_BITS!(int_value >> 26, 0, 6);
        eph_data[2][1] = COMPOSE_BITS!(int_value, 6, 26);
        
        let int_value = Self::unscale_int(ephemeris.w / std::f64::consts::PI, -31);
        eph_data[2][1] |= COMPOSE_BITS!(int_value >> 26, 0, 6);
        eph_data[2][2] = COMPOSE_BITS!(int_value, 6, 26);
        
        let int_value = Self::unscale_int(ephemeris.delta_n / std::f64::consts::PI, -43);
        eph_data[2][2] |= COMPOSE_BITS!(int_value >> 10, 0, 6);
        eph_data[2][3] = COMPOSE_BITS!(int_value, 22, 10);
        
        let int_value = Self::unscale_int(ephemeris.cuc, -29);
        eph_data[2][3] |= COMPOSE_BITS!(int_value, 6, 16);
        
        let int_value = Self::unscale_int(ephemeris.cus, -29);
        eph_data[2][3] |= COMPOSE_BITS!(int_value >> 10, 0, 6);
        eph_data[2][4] = COMPOSE_BITS!(int_value, 22, 10);
        
        let int_value = Self::unscale_int(ephemeris.crc, -5);
        eph_data[2][4] |= COMPOSE_BITS!(int_value, 6, 16);
        
        let int_value = Self::unscale_int(ephemeris.crs, -5);
        eph_data[2][4] |= COMPOSE_BITS!(int_value >> 10, 0, 6);
        eph_data[2][5] = COMPOSE_BITS!(int_value, 22, 10);
        eph_data[2][5] |= COMPOSE_BITS!((ephemeris.toe / 60), 8, 14);
        eph_data[2][6] = 0; // for GST and 8 spare bits

        // PageType 4
        eph_data[3][0] = (4 << 16) | COMPOSE_BITS!(ephemeris.iodc, 6, 10);
        let int_value = Self::unscale_int(ephemeris.cic, -29);
        eph_data[3][0] |= COMPOSE_BITS!(int_value >> 10, 0, 6);
        eph_data[3][1] = COMPOSE_BITS!(int_value, 22, 10);
        
        let int_value = Self::unscale_int(ephemeris.cis, -29);
        eph_data[3][1] |= COMPOSE_BITS!(int_value, 6, 16);
        
        for i in 2..7 {
            eph_data[3][i] = 0; // for GST-UTC, GST-GPS, TOW and 5 spare bits
        }

        0
    }

    fn compose_alm_words(almanac: &[GpsAlmanac], alm_data: &mut [[u32; 7]; 2], week: i32) -> i32 {
        let toa = if !almanac.is_empty() && (almanac[0].valid & 1) != 0 {
            almanac[0].toa
        } else if almanac.len() > 1 && (almanac[1].valid & 1) != 0 {
            almanac[1].toa
        } else if almanac.len() > 2 && (almanac[2].valid & 1) != 0 {
            almanac[2].toa
        } else {
            0
        };

        // Initialize all data to zero
        for i in 0..2 {
            for j in 0..7 {
                alm_data[i][j] = 0;
            }
        }

        if almanac.is_empty() {
            return 0;
        }

        // PageType 5
        alm_data[0][0] = COMPOSE_BITS!(5, 16, 6) | COMPOSE_BITS!(4, 12, 4) | COMPOSE_BITS!(week, 10, 2) | COMPOSE_BITS!((toa / 600), 0, 10);
        alm_data[0][1] = COMPOSE_BITS!(almanac[0].svid, 26, 6); // SVID1 starts here
        
        let int_value = Self::unscale_int(almanac[0].sqrtA - SQRT_A0, -11);
        alm_data[0][1] |= COMPOSE_BITS!(int_value, 13, 13);
        
        let uint_value = Self::unscale_uint(almanac[0].ecc, -16);
        alm_data[0][1] |= COMPOSE_BITS!(uint_value, 2, 11);
        
        let int_value = Self::unscale_int(almanac[0].w / std::f64::consts::PI, -15);
        alm_data[0][1] |= COMPOSE_BITS!(int_value >> 14, 0, 2);
        alm_data[0][2] = COMPOSE_BITS!(int_value, 18, 14);
        
        let int_value = Self::unscale_int((almanac[0].i0 - NOMINAL_I0) / std::f64::consts::PI, -14);
        alm_data[0][2] |= COMPOSE_BITS!(int_value, 7, 11);
        
        let int_value = Self::unscale_int(almanac[0].omega0 / std::f64::consts::PI, -15);
        alm_data[0][2] |= COMPOSE_BITS!(int_value >> 9, 0, 7);
        alm_data[0][3] = COMPOSE_BITS!(int_value, 23, 9);
        
        let int_value = Self::unscale_int(almanac[0].omega_dot / std::f64::consts::PI, -33);
        alm_data[0][3] |= COMPOSE_BITS!(int_value, 12, 11);
        
        let int_value = Self::unscale_int(almanac[0].M0 / std::f64::consts::PI, -15);
        alm_data[0][3] |= COMPOSE_BITS!(int_value >> 4, 0, 12);
        alm_data[0][4] = COMPOSE_BITS!(int_value, 28, 4);
        
        let int_value = Self::unscale_int(almanac[0].af0, -19);
        alm_data[0][4] |= COMPOSE_BITS!(int_value, 12, 16);
        
        let int_value = Self::unscale_int(almanac[0].af1, -38);
        alm_data[0][4] |= COMPOSE_BITS!(int_value >> 1, 0, 12);
        alm_data[0][5] = COMPOSE_BITS!(int_value, 31, 1);

        if almanac.len() > 1 {
            alm_data[0][5] |= COMPOSE_BITS!(if (almanac[1].valid & 1) != 0 { 0 } else { 1 }, 29, 2);
            alm_data[0][5] |= COMPOSE_BITS!(almanac[1].svid, 23, 6); // SVID2 starts here
            
            let int_value = Self::unscale_int(almanac[1].sqrtA - SQRT_A0, -11);
            alm_data[0][5] |= COMPOSE_BITS!(int_value, 10, 13);
            
            let uint_value = Self::unscale_uint(almanac[1].ecc, -16);
            alm_data[0][5] |= COMPOSE_BITS!(uint_value >> 1, 0, 10);
            alm_data[0][6] = COMPOSE_BITS!(uint_value, 31, 1);
            
            let int_value = Self::unscale_int(almanac[1].w / std::f64::consts::PI, -15);
            alm_data[0][6] |= COMPOSE_BITS!(int_value, 15, 16);
            
            let int_value = Self::unscale_int((almanac[1].i0 - NOMINAL_I0) / std::f64::consts::PI, -14);
            alm_data[0][6] |= COMPOSE_BITS!(int_value, 4, 11);
            
            let int_value = Self::unscale_int(almanac[1].omega0 / std::f64::consts::PI, -15);
            alm_data[0][6] |= COMPOSE_BITS!(int_value >> 12, 0, 4);

            // PageType 6
            alm_data[1][0] = COMPOSE_BITS!(6, 16, 6) | COMPOSE_BITS!(4, 12, 4); // Type=6, IODa=4
            alm_data[1][0] |= COMPOSE_BITS!(int_value, 0, 12);
            
            let int_value = Self::unscale_int(almanac[1].omega_dot / std::f64::consts::PI, -33);
            alm_data[1][1] = COMPOSE_BITS!(int_value, 21, 11);
            
            let int_value = Self::unscale_int(almanac[1].M0 / std::f64::consts::PI, -15);
            alm_data[1][1] |= COMPOSE_BITS!(int_value, 5, 16);
            
            let int_value = Self::unscale_int(almanac[1].af0, -19);
            alm_data[1][1] |= COMPOSE_BITS!(int_value >> 11, 0, 5);
            alm_data[1][2] = COMPOSE_BITS!(int_value, 21, 11);
            
            let int_value = Self::unscale_int(almanac[1].af1, -38);
            alm_data[1][2] |= COMPOSE_BITS!(int_value, 8, 13);
            alm_data[1][2] |= COMPOSE_BITS!(if (almanac[1].valid & 1) != 0 { 0 } else { 1 }, 6, 2);

            if almanac.len() > 2 {
                alm_data[1][2] |= COMPOSE_BITS!(almanac[2].svid, 0, 6); // SVID3 starts here
                
                let int_value = Self::unscale_int(almanac[2].sqrtA - SQRT_A0, -11);
                alm_data[1][3] = COMPOSE_BITS!(int_value, 19, 13);
                
                let uint_value = Self::unscale_uint(almanac[2].ecc, -16);
                alm_data[1][3] |= COMPOSE_BITS!(uint_value, 8, 11);
                
                let int_value = Self::unscale_int(almanac[2].w / std::f64::consts::PI, -15);
                alm_data[1][3] |= COMPOSE_BITS!(int_value >> 8, 0, 8);
                alm_data[1][4] = COMPOSE_BITS!(int_value, 24, 8);
                
                let int_value = Self::unscale_int((almanac[2].i0 - NOMINAL_I0) / std::f64::consts::PI, -14);
                alm_data[1][4] |= COMPOSE_BITS!(int_value, 13, 11);
                
                let int_value = Self::unscale_int(almanac[2].omega0 / std::f64::consts::PI, -15);
                alm_data[1][4] |= COMPOSE_BITS!(int_value >> 3, 0, 13);
                alm_data[1][5] = COMPOSE_BITS!(int_value, 29, 3);
                
                let int_value = Self::unscale_int(almanac[2].omega_dot / std::f64::consts::PI, -33);
                alm_data[1][5] |= COMPOSE_BITS!(int_value, 18, 11);
                
                let int_value = Self::unscale_int(almanac[2].M0 / std::f64::consts::PI, -15);
                alm_data[1][5] |= COMPOSE_BITS!(int_value, 2, 16);
                
                let int_value = Self::unscale_int(almanac[2].af0, -19);
                alm_data[1][5] |= COMPOSE_BITS!(int_value >> 14, 0, 2);
                alm_data[1][6] = COMPOSE_BITS!(int_value, 18, 14);
                
                let int_value = Self::unscale_int(almanac[2].af1, -38);
                alm_data[1][6] |= COMPOSE_BITS!(int_value, 5, 13);
                alm_data[1][6] |= COMPOSE_BITS!(if (almanac[2].valid & 1) != 0 { 0 } else { 1 }, 3, 2);
            }
        }

        0
    }

    fn get_page_data(&self, svid: i32, page: i32, subframe: i32, gst: u32, data: &mut [u32; 7]) {
        match page {
            0 => { // page 1
                let tow = gst & 0xFFFFF; // Extract 20-bit TOW
                data.copy_from_slice(&self.gal_eph_data[(svid - 1) as usize][0]);
                // Add TOW (20 bits)
                data[0] |= COMPOSE_BITS!(tow >> 4, 0, 16); // TOW upper 16 bits
                data[1] |= COMPOSE_BITS!(tow, 28, 4); // TOW lower 4 bits
                // add iono correction
                data[3] |= self.gal_iono_data[0];
                data[4] |= self.gal_iono_data[1];
                // add GST (32 bits)
                data[5] |= COMPOSE_BITS!(gst >> 27, 0, 5); // GST upper 5 bits
                data[6] |= COMPOSE_BITS!(gst, 5, 27); // GST lower 27 bits
            },
            1 => { // page 2
                data.copy_from_slice(&self.gal_eph_data[(svid - 1) as usize][1]);
                data[6] = gst;
            },
            2 => { // page 3
                data.copy_from_slice(&self.gal_eph_data[(svid - 1) as usize][2]);
                // add GST
                data[5] |= COMPOSE_BITS!(gst >> 24, 0, 8);
                data[6] |= COMPOSE_BITS!(gst, 8, 24);
            },
            3 => { // page 4
                data.copy_from_slice(&self.gal_eph_data[(svid - 1) as usize][3]);
                // add GST-UTC
                data[1] |= self.gal_utc_data[0];
                data[2] = self.gal_utc_data[1];
                data[3] = self.gal_utc_data[2];
                data[4] = self.gal_utc_data[3];
                // add TOW
                data[6] |= COMPOSE_BITS!(gst, 5, 20);
            },
            4 => { // page 5/6
                data.copy_from_slice(&self.gal_alm_data[(subframe / 2) as usize][(subframe & 1) as usize]);
            },
            _ => {
                // Unknown page, fill with zeros
                for i in 0..7 {
                    data[i] = 0;
                }
            }
        }
    }

    // Helper functions
    fn unscale_int(value: f64, scale: i32) -> i32 {
        let scaled = value * (2.0_f64).powi(scale);
        if scaled >= 0.0 {
            (scaled + 0.5) as i32
        } else {
            (scaled - 0.5) as i32
        }
    }

    fn unscale_uint(value: f64, scale: i32) -> u32 {
        let scaled = value * (2.0_f64).powi(scale);
        (scaled + 0.5) as u32
    }

    fn crc24q_encode(data: &[u32], bit_count: usize) -> u32 {
        crate::crc24q::crc24q_encode(data, bit_count)
    }
}

impl Default for FNavBit {
    fn default() -> Self {
        Self::new()
    }
}