use gnss_rust::json_interpreter::{CNavData, read_nav_file};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== RINEX Parser Test ===");
    println!();
    
    // Путь к тестовому RINEX файлу
    let rinex_file = "/mnt/c/msys64/home/Daniil/gnss_rust/Rinex_Data/rinex_v3_20251560000.rnx";
    
    // Проверяем существование файла
    if !Path::new(rinex_file).exists() {
        eprintln!("ERROR: RINEX file not found: {}", rinex_file);
        return Err("RINEX file not found".into());
    }
    
    println!("[INFO] Testing RINEX file: {}", rinex_file);
    println!();
    
    // Создаем структуру навигационных данных
    let mut nav_data = CNavData::new();
    
    // Читаем RINEX файл
    println!("[INFO] Reading RINEX file...");
    read_nav_file(&mut nav_data, rinex_file);
    println!();
    
    // Выводим статистику
    println!("=== RINEX Parsing Results ===");
    println!("GPS ephemeris count: {}", nav_data.get_gps_ephemeris_count());
    println!("GLONASS ephemeris count: {}", nav_data.get_glonass_ephemeris_count()); 
    println!("BeiDou ephemeris count: {}", nav_data.get_beidou_ephemeris_count());
    println!("Galileo ephemeris count: {}", nav_data.get_galileo_ephemeris_count());
    println!();
    
    // Проверяем ионосферные параметры
    if nav_data.has_gps_iono() {
        println!("[OK] GPS ionospheric parameters found");
    } else {
        println!("[WARNING] GPS ionospheric parameters NOT found");
    }
    
    if nav_data.has_bds_iono() {
        println!("[OK] BeiDou ionospheric parameters found");
    } else {
        println!("[WARNING] BeiDou ionospheric parameters NOT found");
    }
    
    if nav_data.has_gal_iono() {
        println!("[OK] Galileo ionospheric parameters found");
    } else {
        println!("[WARNING] Galileo ionospheric parameters NOT found");
    }
    println!();
    
    // Проверяем UTC параметры
    if nav_data.has_gps_utc() {
        println!("[OK] GPS UTC parameters found");
    } else {
        println!("[WARNING] GPS UTC parameters NOT found");
    }
    println!();
    
    // Показываем первые несколько эфемерид для проверки
    if nav_data.get_gps_ephemeris_count() > 0 {
        println!("=== First GPS Ephemeris Sample ===");
        if let Some(first_eph) = nav_data.get_first_gps_ephemeris() {
            println!("SVID: {}", first_eph.svid);
            println!("TOE: {}", first_eph.toe);
            println!("TOC: {}", first_eph.toc);
            println!("sqrt(A): {:.6}", first_eph.sqrtA);
            println!("Eccentricity: {:.9}", first_eph.ecc);
            println!("M0: {:.9}", first_eph.M0);
            println!("af0: {:.12e}", first_eph.af0);
            println!("af1: {:.12e}", first_eph.af1);
            println!("af2: {:.12e}", first_eph.af2);
        }
        println!();
    }
    
    println!("=== Test Completed ===");
    Ok(())
}