use gnss_rust::json_interpreter::{CNavData, read_nav_file_limited};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== RINEX Data Validation Test ===");
    println!();
    
    let rinex_file = "/mnt/c/msys64/home/Daniil/gnss_rust/Rinex_Data/rinex_v3_20251560000.rnx";
    
    if !Path::new(rinex_file).exists() {
        eprintln!("ERROR: RINEX file not found: {}", rinex_file);
        return Err("RINEX file not found".into());
    }
    
    println!("[INFO] Validating RINEX data structures in: {}", rinex_file);
    println!("[INFO] Reading limited sample from each GNSS system for validation");
    println!();
    
    let mut nav_data = CNavData::new();
    read_nav_file_limited(&mut nav_data, rinex_file, 3); // Читаем по 3 эфемериды от каждой системы
    
    // Проверяем GPS эфемериды
    if let Some(first_gps) = nav_data.get_first_gps_ephemeris() {
        println!("=== GPS Ephemeris Validation (First Entry) ===");
        println!("SVID: {}", first_gps.svid);
        println!("Week: {}", first_gps.week);
        println!("TOE: {} sec", first_gps.toe);
        println!("TOC: {} sec", first_gps.toc);
        println!("sqrt(A): {:.6} m^1/2", first_gps.sqrtA);
        println!("Semi-major axis: {:.3} km", first_gps.axis / 1000.0);
        println!("Eccentricity: {:.9}", first_gps.ecc);
        println!("Inclination (i0): {:.9} rad ({:.3}°)", first_gps.i0, first_gps.i0.to_degrees());
        println!("RAAN (Omega0): {:.9} rad ({:.3}°)", first_gps.omega0, first_gps.omega0.to_degrees());
        println!("Arg of Perigee (w): {:.9} rad ({:.3}°)", first_gps.w, first_gps.w.to_degrees());
        println!("Mean Anomaly (M0): {:.9} rad ({:.3}°)", first_gps.M0, first_gps.M0.to_degrees());
        println!("Mean Motion Diff (Δn): {:.12} rad/s", first_gps.delta_n);
        println!("RAAN Rate (Ω̇): {:.12} rad/s", first_gps.omega_dot);
        println!("Inclination Rate (İ): {:.12} rad/s", first_gps.idot);
        println!();
        
        // Корректировочные коэффициенты
        println!("=== Correction Coefficients ===");
        println!("Crs: {:.6} m", first_gps.crs);
        println!("Crc: {:.6} m", first_gps.crc);
        println!("Cus: {:.9} rad", first_gps.cus);
        println!("Cuc: {:.9} rad", first_gps.cuc);
        println!("Cis: {:.9} rad", first_gps.cis);
        println!("Cic: {:.9} rad", first_gps.cic);
        println!();
        
        // Часовые поправки
        println!("=== Clock Corrections ===");
        println!("af0: {:.12e} s", first_gps.af0);
        println!("af1: {:.12e} s/s", first_gps.af1);
        println!("af2: {:.12e} s/s²", first_gps.af2);
        println!("TGD: {:.12e} s", first_gps.tgd);
        println!();
        
        // Статусные параметры
        println!("=== Status Parameters ===");
        println!("IODE: {}", first_gps.iode);
        println!("IODC: {}", first_gps.iodc);
        println!("URA: {}", first_gps.ura);
        println!("Health: {}", first_gps.health);
        println!("Valid: {}", first_gps.valid);
        println!("Flag: {}", first_gps.flag);
        println!();
        
        // Проверим физические ограничения
        let mut errors = Vec::new();
        
        // Проверка орбитальных параметров
        if first_gps.sqrtA < 5000.0 || first_gps.sqrtA > 6500.0 {
            errors.push(format!("Неправильная sqrt(A): {:.1} (должна быть 5000-6500)", first_gps.sqrtA));
        }
        
        if first_gps.ecc < 0.0 || first_gps.ecc > 0.1 {
            errors.push(format!("Неправильная эксцентричность: {:.6} (должна быть 0-0.1)", first_gps.ecc));
        }
        
        let inclination_deg = first_gps.i0.to_degrees();
        if inclination_deg < 50.0 || inclination_deg > 60.0 {
            errors.push(format!("Неправильная наклонность: {:.1}° (должна быть 50-60°)", inclination_deg));
        }
        
        // Проверка часовых поправок
        if first_gps.af0.abs() > 1e-3 {
            errors.push(format!("Слишком большая af0: {:.3e} (должна быть < 1e-3)", first_gps.af0));
        }
        
        if errors.is_empty() {
            println!("✅ GPS ephemeris data looks physically reasonable");
        } else {
            println!("❌ GPS ephemeris validation errors:");
            for error in &errors {
                println!("   - {}", error);
            }
        }
    } else {
        println!("❌ No GPS ephemeris found in parsed data");
    }
    
    println!();
    println!("{}", "=".repeat(80));
    
    // Проверяем GLONASS эфемериды
    if nav_data.get_glonass_ephemeris_count() > 0 {
        println!("=== GLONASS Ephemeris Sample ===");
        let glonass_list = nav_data.get_glonass_ephemeris();
        if let Some(first_glo) = glonass_list.first() {
            println!("Slot: {}", first_glo.slot);
            println!("Frequency number: {}", first_glo.freq);
            println!("tk (time mark): {} ({}h {}m)", 
                     first_glo.tk,
                     (first_glo.tk >> 7) & 0x1F,
                     (first_glo.tk >> 1) & 0x3F);
            println!("Position: [{:.3}, {:.3}, {:.3}] km", 
                     first_glo.x / 1000.0, first_glo.y / 1000.0, first_glo.z / 1000.0);
            println!("Velocity: [{:.3}, {:.3}, {:.3}] km/s", 
                     first_glo.vx / 1000.0, first_glo.vy / 1000.0, first_glo.vz / 1000.0);
            println!("Clock bias (tn): {:.9e} s", -first_glo.tn);
            println!("Relative freq bias (gamma): {:.9e}", first_glo.gamma);
            println!("✅ GLONASS ephemeris data structure looks good");
        }
    } else {
        println!("❌ No GLONASS ephemeris found in parsed data");
    }
    
    println!();
    println!("{}", "=".repeat(80));
    
    // Проверяем BeiDou эфемериды
    if nav_data.get_beidou_ephemeris_count() > 0 {
        println!("=== BeiDou Ephemeris Sample ===");
        let beidou_list = nav_data.get_beidou_ephemeris();
        if let Some(first_bds) = beidou_list.first() {
            println!("SVID: {}", first_bds.svid);
            println!("TOE: {} sec (BDT corrected)", first_bds.toe);
            println!("sqrt(A): {:.6} m^1/2", first_bds.sqrtA);
            println!("Eccentricity: {:.9}", first_bds.ecc);
            println!("Inclination: {:.3}°", first_bds.i0.to_degrees());
            println!("✅ BeiDou ephemeris data structure looks good");
        }
    } else {
        println!("❌ No BeiDou ephemeris found in parsed data");
    }
    
    println!();
    println!("{}", "=".repeat(80));
    
    // Проверяем Galileo эфемериды  
    if nav_data.get_galileo_ephemeris_count() > 0 {
        println!("=== Galileo Ephemeris Sample ===");
        let galileo_list = nav_data.get_galileo_ephemeris();
        if let Some(first_gal) = galileo_list.first() {
            println!("SVID: {}", first_gal.svid);
            println!("TOE: {} sec", first_gal.toe);
            println!("sqrt(A): {:.6} m^1/2", first_gal.sqrtA);
            println!("Eccentricity: {:.9}", first_gal.ecc);
            println!("Inclination: {:.3}°", first_gal.i0.to_degrees());
            println!("✅ Galileo ephemeris data structure looks good");
        }
    } else {
        println!("❌ No Galileo ephemeris found in parsed data");
    }
    
    println!();
    println!("=== Summary ===");
    println!("Total GPS ephemeris loaded: {}", nav_data.get_gps_ephemeris_count());
    println!("Total GLONASS ephemeris loaded: {}", nav_data.get_glonass_ephemeris_count());
    println!("Total BeiDou ephemeris loaded: {}", nav_data.get_beidou_ephemeris_count());
    println!("Total Galileo ephemeris loaded: {}", nav_data.get_galileo_ephemeris_count());
    
    // Проверяем ионосферные параметры
    if nav_data.has_gps_iono() {
        println!("✅ GPS ionospheric parameters are present");
    } else {
        println!("❌ GPS ionospheric parameters are missing");
    }
    
    if nav_data.has_bds_iono() {
        println!("✅ BeiDou ionospheric parameters are present");
    } else {
        println!("❌ BeiDou ionospheric parameters are missing");
    }
    
    if nav_data.has_gal_iono() {
        println!("✅ Galileo ionospheric parameters are present");
    } else {
        println!("❌ Galileo ionospheric parameters are missing");
    }
    
    println!();
    println!("=== Validation Complete ===");
    
    Ok(())
}