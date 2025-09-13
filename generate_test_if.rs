use std::fs::File;
use std::io::Write;

fn main() -> std::io::Result<()> {
    println!("Генерация тестового IF файла с GNSS-подобными сигналами...");
    
    let sample_rate = 5_000_000.0; // 5 МГц
    let duration_ms = 10; // 10 миллисекунд
    let samples_per_ms = (sample_rate / 1000.0) as usize;
    let total_samples = samples_per_ms * duration_ms;
    
    // Симуляция 5 спутниковых сигналов с разными доплеровскими сдвигами
    let satellites = vec![
        (3500.0, 0.1),   // Спутник 1: +3.5 кГц доплер, амплитуда 0.1
        (-2100.0, 0.08), // Спутник 2: -2.1 кГц доплер
        (1200.0, 0.12),  // Спутник 3: +1.2 кГц доплер
        (-4800.0, 0.07), // Спутник 4: -4.8 кГц доплер
        (500.0, 0.09),   // Спутник 5: +0.5 кГц доплер
    ];
    
    let mut iq8_data = Vec::with_capacity(total_samples * 2);
    
    for sample_idx in 0..total_samples {
        let t = sample_idx as f64 / sample_rate;
        let mut i_total = 0.0;
        let mut q_total = 0.0;
        
        // Суммируем сигналы всех спутников
        for &(doppler_hz, amplitude) in &satellites {
            let phase = 2.0 * std::f64::consts::PI * doppler_hz * t;
            i_total += amplitude * phase.cos();
            q_total += amplitude * phase.sin();
        }
        
        // Добавляем шум
        i_total += (rand::random::<f64>() - 0.5) * 0.02;
        q_total += (rand::random::<f64>() - 0.5) * 0.02;
        
        // Квантуем в IQ8 формат (8 бит на компоненту)
        let i_quant = (i_total * 127.0).max(-128.0).min(127.0) as i8;
        let q_quant = (q_total * 127.0).max(-128.0).min(127.0) as i8;
        
        iq8_data.push(i_quant as u8);
        iq8_data.push(q_quant as u8);
    }
    
    // Сохраняем в файл
    let mut file = File::create("test_if_data.bin")?;
    file.write_all(&iq8_data)?;
    
    println!("Создан файл test_if_data.bin:");
    println!("  - Формат: IQ8");
    println!("  - Частота дискретизации: {} МГц", sample_rate / 1e6);
    println!("  - Длительность: {} мс", duration_ms);
    println!("  - Размер: {} байт", iq8_data.len());
    println!("  - Спутников: {}", satellites.len());
    
    Ok(())
}

// Простая реализация rand для теста
mod rand {
    pub fn random<T>() -> T 
    where T: Default {
        // Простой псевдослучайный генератор
        use std::time::{SystemTime, UNIX_EPOCH};
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .subsec_nanos();
        
        // Для f64
        if std::any::TypeId::of::<T>() == std::any::TypeId::of::<f64>() {
            let value = (nanos as f64) / (u32::MAX as f64);
            unsafe { std::mem::transmute_copy(&value) }
        } else {
            T::default()
        }
    }
}