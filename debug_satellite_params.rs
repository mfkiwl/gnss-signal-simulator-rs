// Тест для проверки заполнения параметров спутника

use gnss_rust::*;
use gnss_rust::types::*;
use gnss_rust::json_interpreter::*;
use gnss_rust::satellite_param::*;
use gnss_rust::gnsstime::*;
use gnss_rust::coordinate::*;
use std::path::Path;

fn main() {
    println!("🔍 Проверка заполнения параметров спутника");
    println!("===========================================");
    
    // Загружаем конфигурацию
    let preset_path = "presets/GPS_BDS_GAL_triple_system.json";
    if !Path::new(preset_path).exists() {
        eprintln!("❌ Файл пресета не найден: {}", preset_path);
        return;
    }
    
    // Создаем экземпляр генератора и парсим настройки
    let mut generator = json_interpreter::IFDataGenSettings::new();
    match generator.load_from_json(preset_path) {
        Ok(nav_data) => {
            println!("✅ Конфигурация загружена успешно");
            
            // Тестовые данные приемника (Chicago)
            let receiver_pos = KinematicInfo {
                x: -3.977830e6,
                y: -5.165207e6, 
                z: 4.221499e6,
                vx: 0.0, vy: 0.0, vz: 0.0,
            };
            
            let receiver_lla = ecef_to_lla(&receiver_pos);
            println!("📍 Позиция приемника: {:.3}° N, {:.3}° W, {:.1} м", 
                     receiver_lla.lat.to_degrees(), 
                     -receiver_lla.lon.to_degrees(), 
                     receiver_lla.alt);
            
            // Тестовое время
            let test_time = GnssTime {
                Week: 2369,
                MilliSeconds: 36330000, // 10:05:30
                SubMilliSeconds: 0.0,
            };
            
            // Тестируем GPS спутник
            if let Some(gps_ephs) = nav_data.get_gps_ephemeris() {
                if !gps_ephs.is_empty() {
                    let eph = &gps_ephs[0];
                    println!("\n🛰️  Тестирование GPS G{:02}:", eph.svid);
                    
                    let mut sat_param = SatelliteParam::default();
                    let iono_param = nav_data.get_gps_iono().cloned().unwrap_or_default();
                    
                    // Вызываем функцию расчета параметров
                    get_satellite_param(
                        &receiver_pos,
                        &receiver_lla,
                        &test_time,
                        GnssSystem::GpsSystem,
                        eph,
                        &iono_param,
                        &mut sat_param
                    );
                    
                    // Проверяем заполненность параметров
                    println!("   SVID: {}", sat_param.svid);
                    println!("   Система: {:?}", sat_param.system);
                    println!("   Позиция: ({:.1}, {:.1}, {:.1}) км", 
                             sat_param.PosVel.x/1000.0, sat_param.PosVel.y/1000.0, sat_param.PosVel.z/1000.0);
                    println!("   Скорость: ({:.3}, {:.3}, {:.3}) км/с", 
                             sat_param.PosVel.vx/1000.0, sat_param.PosVel.vy/1000.0, sat_param.PosVel.vz/1000.0);
                    println!("   Угол места: {:.1}°", sat_param.Elevation.to_degrees());
                    println!("   Азимут: {:.1}°", sat_param.Azimuth.to_degrees());
                    println!("   Радиальная скорость: {:.1} м/с", sat_param.RelativeSpeed);
                    println!("   Время распространения: {:.6} с", sat_param.TravelTime);
                    println!("   Ионосферная задержка: {:.6} м", sat_param.IonoDelay);
                    println!("   CN0: {} дБ-Гц", sat_param.CN0 / 100);
                    println!("   Временная метка: {}", sat_param.PosTimeTag);
                    
                    // Проверяем доплеровский сдвиг
                    let doppler = get_doppler(&sat_param, SIGNAL_INDEX_L1CA);
                    println!("   🎯 Доплеровский сдвиг L1CA: {:.1} Гц", doppler);
                    
                    // Проверки корректности
                    let mut checks_passed = 0;
                    let total_checks = 8;
                    
                    if sat_param.svid == eph.svid as i32 { checks_passed += 1; }
                    else { println!("   ❌ SVID не совпадает"); }
                    
                    if sat_param.PosVel.x.abs() > 1e6 { checks_passed += 1; }
                    else { println!("   ❌ Позиция X не заполнена"); }
                    
                    if sat_param.PosVel.y.abs() > 1e6 { checks_passed += 1; }
                    else { println!("   ❌ Позиция Y не заполнена"); }
                    
                    if sat_param.PosVel.z.abs() > 1e6 { checks_passed += 1; }
                    else { println!("   ❌ Позиция Z не заполнена"); }
                    
                    if sat_param.Elevation > 0.0 { checks_passed += 1; }
                    else { println!("   ❌ Угол места не заполнен или отрицательный"); }
                    
                    if sat_param.RelativeSpeed.abs() > 1.0 { checks_passed += 1; }
                    else { println!("   ❌ Радиальная скорость слишком мала"); }
                    
                    if sat_param.TravelTime > 0.0 { checks_passed += 1; }
                    else { println!("   ❌ Время распространения не заполнено"); }
                    
                    if sat_param.PosTimeTag != 0 { checks_passed += 1; }
                    else { println!("   ❌ Временная метка не заполнена"); }
                    
                    println!("   📊 Проверок пройдено: {}/{}", checks_passed, total_checks);
                    
                    if checks_passed == total_checks {
                        println!("   ✅ ВСЕ ПАРАМЕТРЫ ЗАПОЛНЕНЫ КОРРЕКТНО!");
                    } else {
                        println!("   ⚠️  Обнаружены проблемы с заполнением параметров");
                    }
                }
            }
        },
        Err(e) => {
            eprintln!("❌ Ошибка загрузки конфигурации: {:?}", e);
        }
    }
}