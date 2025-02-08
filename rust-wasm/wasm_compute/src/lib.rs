use wasm_bindgen::prelude::*;

/// ALGORITHM 1: Integer Matrix Multiplication  
/// (Multiply two 400×400 matrices with modulo arithmetic)
#[wasm_bindgen]
pub fn algo1() -> u32 {
    let size: usize = 400;
    let mod_val: u32 = 10007;
    let mut matrix_a = vec![0u32; size * size];
    let mut matrix_b = vec![0u32; size * size];
    let mut result = vec![0u32; size * size];
    
    for i in 0..size {
        for j in 0..size {
            let i_u32 = i as u32;
            let j_u32 = j as u32;   
            matrix_a[i * size + j] = ((i_u32 + j_u32) % mod_val) as u32;
            matrix_b[i * size + j] = ((i_u32 * j_u32) % mod_val) as u32;
        }
    }
    
    for i in 0..size {
        for j in 0..size {
            let mut sum = 0u32;
            for k in 0..size {
                sum = (sum + matrix_a[i * size + k] * matrix_b[k * size + j]) % mod_val;
            }
            result[i * size + j] = sum;
        }
    }
    
    result.iter().fold(0, |acc, &x| (acc + x) % mod_val)
}

/// ALGORITHM 2: Fibonacci Iterative  
/// (Loop 10,000,000 times computing Fibonacci modulo 1,000,000,007)
#[wasm_bindgen]
pub fn algo2() -> u32 {
    let iterations = 10_000_000;
    let mod_val = 1_000_000_007;
    let mut a: u32 = 0;
    let mut b: u32 = 1;
    for _ in 0..iterations {
        let c = (a + b) % mod_val;
        a = b;
        b = c;
    }
    a
}

/// ALGORITHM 3: Sieve of Eratosthenes  
/// (Compute primes up to 50,000 and sum them modulo 10007)
#[wasm_bindgen]
pub fn algo3() -> u32 {
    let limit = 50000;
    let mod_val = 10007;
    let mut is_prime = vec![true; (limit + 1) as usize];
    if limit >= 0 { is_prime[0] = false; }
    if limit >= 1 { is_prime[1] = false; }
    let r = (limit as f64).sqrt() as u32;
    for i in 2..=r {
        if is_prime[i as usize] {
            let mut j = i * i;
            while j <= limit {
                is_prime[j as usize] = false;
                j += i;
            }
        }
    }
    let mut sum = 0;
    for i in 2..=limit {
        if is_prime[i as usize] {
            sum = (sum + i) % mod_val;
        }
    }
    sum
}

/// ALGORITHM 4: Sorting a Large Array  
/// (Generate 100,000 pseudo‑random numbers, sort them, and sum modulo 10007)
#[wasm_bindgen]
pub fn algo4() -> u32 {
    let size = 100_000;
    let mod_val = 10007;
    let mut arr: Vec<u32> = Vec::with_capacity(size);
    let mut seed = 42u32;
    for _ in 0..size {
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        arr.push(seed % 100000);
    }
    arr.sort_unstable();
    arr.iter().fold(0, |acc, &x| (acc + x) % mod_val)
}

/// ALGORITHM 5: Mandelbrot Set Calculation  
/// (For a 300×300 grid, count points that remain bounded after 100 iterations)
#[wasm_bindgen]
pub fn algo5() -> u32 {
    let width = 300;
    let height = 300;
    let max_iter = 100;
    let mut count = 0;
    for i in 0..width {
        for j in 0..height {
            let cx = (i as f64 / width as f64) * 3.5 - 2.5;
            let cy = (j as f64 / height as f64) * 2.0 - 1.0;
            let (mut x, mut y) = (0.0, 0.0);
            let mut iter = 0;
            while x*x + y*y <= 4.0 && iter < max_iter {
                let temp = x*x - y*y + cx;
                y = 2.0*x*y + cy;
                x = temp;
                iter += 1;
            }
            if iter == max_iter {
                count += 1;
            }
        }
    }
    count
}

/// ALGORITHM 6: Monte Carlo Pi Estimation  
/// (Estimate π using 1,000,000 pseudo‑random points)
#[wasm_bindgen]
pub fn algo6() -> u32 {
    let iterations = 1_000_000;
    let mut inside = 0;
    let mut seed = 42u32;
    for _ in 0..iterations {
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        let x = (seed % 10000) as f64 / 10000.0;
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        let y = (seed % 10000) as f64 / 10000.0;
        if x*x + y*y <= 1.0 {
            inside += 1;
        }
    }
    let pi_est = 4.0 * inside as f64 / iterations as f64;
    (pi_est * 1000.0) as u32
}

/// ALGORITHM 7: N‑Queens Solver (N = 10)  
/// (Count the number of solutions using backtracking)
fn solve_n_queens(n: usize) -> u32 {
    fn solve(row: usize, n: usize, cols: &mut Vec<bool>, diag1: &mut Vec<bool>, diag2: &mut Vec<bool>, count: &mut u32) {
        if row == n {
            *count += 1;
            return;
        }
        for col in 0..n {
            if !cols[col] && !diag1[row+col] && !diag2[row+n-col-1] {
                cols[col] = true;
                diag1[row+col] = true;
                diag2[row+n-col-1] = true;
                solve(row+1, n, cols, diag1, diag2, count);
                cols[col] = false;
                diag1[row+col] = false;
                diag2[row+n-col-1] = false;
            }
        }
    }
    let mut count = 0;
    let mut cols = vec![false; n];
    let mut diag1 = vec![false; 2*n - 1];
    let mut diag2 = vec![false; 2*n - 1];
    solve(0, n, &mut cols, &mut diag1, &mut diag2, &mut count);
    count
}

#[wasm_bindgen]
pub fn algo7() -> u32 {
    solve_n_queens(10)
}

/// ALGORITHM 8: Factorial Modulo  
/// (Compute 20000! modulo 1,000,000,007)
#[wasm_bindgen]
pub fn algo8() -> u32 {
    let n = 20000;
    let mod_val = 1_000_000_007u32;
    let mut fact = 1u32;
    for i in 1..=n {
        fact = fact.wrapping_mul(i) % mod_val;
    }
    fact
}

/// ALGORITHM 9: GCD Sum  
/// (For i = 1 to 100000, compute gcd(i, i+100) and sum modulo 10007)
fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let temp = a % b;
        a = b;
        b = temp;
    }
    a
}

#[wasm_bindgen]
pub fn algo9() -> u32 {
    let limit = 100000;
    let mod_val = 10007;
    let mut sum = 0u32;
    for i in 1..=limit {
        let g = gcd(i, i + 100);
        sum = (sum + g) % mod_val;
    }
    sum
}

/// ALGORITHM 10: Discrete Fourier Transform (DFT)  
/// (Compute a naive DFT on a vector of size 256 and sum the magnitudes)
#[wasm_bindgen]
pub fn algo10() -> u32 {
    let n = 256;
    let input: Vec<f64> = (0..n).map(|i| (i as f64).sin()).collect();
    let mut sum = 0f64;
    for k in 0..n {
        let (mut re, mut im) = (0f64, 0f64);
        for t in 0..n {
            let angle = 2.0 * std::f64::consts::PI * (t * k) as f64 / n as f64;
            re += input[t] * angle.cos();
            im -= input[t] * angle.sin();
        }
        sum += (re*re + im*im).sqrt();
    }
    sum as u32
}

/// ALGORITHM 11: 1D Convolution  
/// (Convolve an array of 1000 elements with a kernel of size 50 and sum modulo 10007)
#[wasm_bindgen]
pub fn algo11() -> u32 {
    let size = 1000;
    let kernel_size = 50;
    let mod_val = 10007;
    let mut data = vec![0u32; size];
    let mut kernel = vec![0u32; kernel_size];
    for i in 0..size {
        data[i] = (i % 100) as u32;
    }
    for i in 0..kernel_size {
        kernel[i] = ((i+1) % 10) as u32;
    }
    let mut result = vec![0u32; size];
    for i in 0..size {
        let mut sum = 0u32;
        for j in 0..kernel_size {
            if i >= j {
                sum += data[i - j] * kernel[j];
            }
        }
        result[i] = sum % mod_val;
    }
    result.iter().fold(0, |acc, &x| (acc + x) % mod_val)
}

/// ALGORITHM 12: Dot Product  
/// (Compute the dot product of two arrays of 1,000,000 elements modulo 10007)
#[wasm_bindgen]
pub fn algo12() -> u32 {
    let size = 1_000_000;
    let mod_val = 10007;
    let vec1: Vec<u32> = (0..size).map(|i| (i % 1000) as u32).collect();
    let vec2: Vec<u32> = (0..size).map(|i| ((size - i) % 1000) as u32).collect();
    let mut dot = 0u32;
    for i in 0..size {
        dot = (dot + vec1[i] * vec2[i]) % mod_val;
    }
    dot
}

/// ALGORITHM 13: Harmonic Sum  
/// (Compute the harmonic series sum up to 10,000,000 and multiply by 1000)
#[wasm_bindgen]
pub fn algo13() -> u32 {
    let n = 10_000_000;
    let mut sum = 0f64;
    for i in 1..=n {
        sum += 1.0 / (i as f64);
    }
    (sum * 1000.0) as u32
}

/// ALGORITHM 14: Bubble Sort  
/// (Sort an array of 1000 numbers using bubble sort and sum modulo 10007)
#[wasm_bindgen]
pub fn algo14() -> u32 {
    let size = 1000;
    let mod_val = 10007;
    let mut arr: Vec<u32> = (0..size).map(|i| (size - i) as u32).collect();
    for _ in 0..size {
        for j in 0..size-1 {
            if arr[j] > arr[j+1] {
                arr.swap(j, j+1);
            }
        }
    }
    arr.iter().fold(0, |acc, &x| (acc + x) % mod_val)
}

/// ALGORITHM 15: Maximum Subarray Sum (Kadane's Algorithm)  
/// (Compute the maximum subarray sum on an array of 100,000 pseudo‑random integers)
#[wasm_bindgen]
pub fn algo15() -> u32 {
    let size = 100000;
    let mut arr = vec![0i32; size];
    let mut seed = 42u32;
    for i in 0..size {
        seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
        arr[i] = (seed % 1001) as i32 - 500;
    }
    let mut max_current = arr[0];
    let mut max_global = arr[0];
    for i in 1..size {
        max_current = std::cmp::max(arr[i], max_current + arr[i]);
        if max_current > max_global {
            max_global = max_current;
        }
    }
    if max_global < 0 { 0 } else { max_global as u32 }
}
