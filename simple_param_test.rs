// Простая проверка заполнения параметров спутника

use gnss_rust::types::*;
use gnss_rust::satellite_param::*;
use gnss_rust::coordinate::*;

fn main() {
    println!("🔍 Тест заполнения параметров SatelliteParam");
    println!("============================================");
    
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
    
    // Тестовые эфемериды GPS (упрощенные)
    let test_eph = GpsEphemeris {
        svid: 1,
        week: 2369,
        toc: 36330.0,
        toe: 36330.0,
        M0: 2.0,
        delta_n: 0.0,
        ecc: 0.01,
        sqrtA: 5153.0, // типичное значение для GPS
        omega0: 1.5,
        i0: 0.97, // ~55 градусов
        omega: 0.5,
        omega_dot: -2.6e-9,
        idot: 0.0,
        cuc: 0.0,
        cus: 0.0,
        crc: 0.0,
        crs: 0.0,
        cic: 0.0,
        cis: 0.0,
        af0: 0.0,
        af1: 0.0,
        af2: 0.0,
        tgd: 0.0,
        iodc: 100,
        iode: 100,
        health: 0,
        ura: 0,
        code_on_l2: 0,
        l2_p_data_flag: 0,
        fit_interval: 4,
        aodo: 0,
        omega_t: 1.5,
        omega_delta: -2.6e-9,
        Ek: 0.0,
        Ek_dot: 0.0,
        tgd_ext: [0.0; 5],
        tgd2: 0.0,
    };
    
    // Тестовые ионосферные параметры
    let iono_param = IonoParam {
        a0: 1.4e-8,
        a1: 1.4e-8, 
        a2: -6.0e-8,
        a3: -1.2e-7,
        b0: 1.0e5,
        b1: 9.8e4,
        b2: -6.5e4,
        b3: -3.3e5,
        flag: 0,
    };
    
    // Создаем структуру параметров спутника
    let mut sat_param = SatelliteParam::default();
    
    println!("\n🛰️  Тестирование расчета параметров для GPS G01:");
    println!("📊 ПЕРЕД вызовом get_satellite_param:");
    println!("   PosVel.x: {:.3} м", sat_param.PosVel.x);
    println!("   PosVel.y: {:.3} м", sat_param.PosVel.y);
    println!("   PosVel.z: {:.3} м", sat_param.PosVel.z);
    println!("   PosTimeTag: {}", sat_param.PosTimeTag);
    println!("   RelativeSpeed: {:.3} м/с", sat_param.RelativeSpeed);
    println!("   Elevation: {:.3}°", sat_param.Elevation.to_degrees());
    println!("   Azimuth: {:.3}°", sat_param.Azimuth.to_degrees());
    println!("   TravelTime: {:.6} с", sat_param.TravelTime);
    
    // Вызываем функцию расчета параметров
    get_satellite_param(
        &receiver_pos,
        &receiver_lla,
        &test_time,
        GnssSystem::GpsSystem,
        &test_eph,
        &iono_param,
        &mut sat_param
    );
    
    println!("\n📊 ПОСЛЕ вызова get_satellite_param:");
    println!("   SVID: {}", sat_param.svid);
    println!("   Система: {:?}", sat_param.system);
    println!("   PosVel.x: {:.3} м", sat_param.PosVel.x);
    println!("   PosVel.y: {:.3} м", sat_param.PosVel.y);
    println!("   PosVel.z: {:.3} м", sat_param.PosVel.z);
    println!("   PosVel.vx: {:.3} м/с", sat_param.PosVel.vx);
    println!("   PosVel.vy: {:.3} м/с", sat_param.PosVel.vy);
    println!("   PosVel.vz: {:.3} м/с", sat_param.PosVel.vz);
    println!("   PosTimeTag: {}", sat_param.PosTimeTag);
    println!("   Elevation: {:.3}°", sat_param.Elevation.to_degrees());
    println!("   Azimuth: {:.3}°", sat_param.Azimuth.to_degrees());
    println!("   RelativeSpeed: {:.3} м/с", sat_param.RelativeSpeed);
    println!("   TravelTime: {:.6} с", sat_param.TravelTime);
    println!("   IonoDelay: {:.6} м", sat_param.IonoDelay);
    println!("   CN0: {} дБ-Гц", sat_param.CN0);
    
    // Проверяем доплеровский сдвиг
    let doppler = get_doppler(&sat_param, gnss_rust::SIGNAL_INDEX_L1CA);
    println!("   🎯 Доплеровский сдвиг L1CA: {:.1} Гц", doppler);
    
    // Проверки корректности
    println!("\n✅ Проверки корректности:");
    let mut checks_passed = 0;
    let total_checks = 9;
    
    if sat_param.svid == test_eph.svid as i32 { 
        println!("   ✓ SVID корректен: {}", sat_param.svid);
        checks_passed += 1; 
    } else { 
        println!("   ❌ SVID не совпадает: {} != {}", sat_param.svid, test_eph.svid); 
    }
    
    if sat_param.PosVel.x.abs() > 1e6 { 
        println!("   ✓ Позиция X заполнена: {:.1} км", sat_param.PosVel.x/1000.0);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Позиция X не заполнена: {}", sat_param.PosVel.x); 
    }
    
    if sat_param.PosVel.y.abs() > 1e6 { 
        println!("   ✓ Позиция Y заполнена: {:.1} км", sat_param.PosVel.y/1000.0);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Позиция Y не заполнена: {}", sat_param.PosVel.y); 
    }
    
    if sat_param.PosVel.z.abs() > 1e6 { 
        println!("   ✓ Позиция Z заполнена: {:.1} км", sat_param.PosVel.z/1000.0);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Позиция Z не заполнена: {}", sat_param.PosVel.z); 
    }
    
    if sat_param.Elevation > 0.0 { 
        println!("   ✓ Угол места положительный: {:.1}°", sat_param.Elevation.to_degrees());
        checks_passed += 1; 
    } else { 
        println!("   ❌ Угол места не заполнен или отрицательный: {:.1}°", sat_param.Elevation.to_degrees()); 
    }
    
    if sat_param.RelativeSpeed.abs() > 1.0 { 
        println!("   ✓ Радиальная скорость разумная: {:.1} м/с", sat_param.RelativeSpeed);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Радиальная скорость слишком мала: {:.1} м/с", sat_param.RelativeSpeed); 
    }
    
    if sat_param.TravelTime > 0.0 { 
        println!("   ✓ Время распространения положительное: {:.6} с", sat_param.TravelTime);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Время распространения не заполнено: {:.6} с", sat_param.TravelTime); 
    }
    
    if sat_param.PosTimeTag != 0 { 
        println!("   ✓ Временная метка установлена: {}", sat_param.PosTimeTag);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Временная метка не заполнена: {}", sat_param.PosTimeTag); 
    }
    
    if doppler.abs() > 0.1 { 
        println!("   ✓ Доплеровский сдвиг разумный: {:.1} Гц", doppler);
        checks_passed += 1; 
    } else { 
        println!("   ❌ Доплеровский сдвиг слишком мал: {:.1} Гц", doppler); 
    }
    
    println!("\n📊 ИТОГИ:");
    println!("   Проверок пройдено: {}/{}", checks_passed, total_checks);
    
    if checks_passed == total_checks {
        println!("   🎉 ВСЕ ПАРАМЕТРЫ ЗАПОЛНЕНЫ КОРРЕКТНО!");
    } else {
        println!("   ⚠️  Обнаружены проблемы с заполнением {} параметров", total_checks - checks_passed);
    }
}