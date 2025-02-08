// src/App.tsx
import { useState } from 'react';
import init, {
  algo1, algo2, algo3, algo4, algo5,
  algo6, algo7, algo8, algo9, algo10,
  algo11, algo12, algo13, algo14, algo15,
} from '../rust-wasm/wasm_compute/pkg/wasm_compute';

function App() {
  // Create an array of 15 objects to store timings and results.
  const initialResults = Array.from({ length: 15 }, () => ({
    wasmTime: null as number | null,
    wasmResult: null as number | null,
    jsTime: null as number | null,
    jsResult: null as number | null,
  }));
  const [results, setResults] = useState(initialResults);

  // Labels for each algorithm.
  const algoLabels = [
    "Matrix Multiplication",
    "Fibonacci Iterative",
    "Sieve of Eratosthenes",
    "Sort Array",
    "Mandelbrot Set",
    "Monte Carlo Pi",
    "N-Queens",
    "Factorial Modulo",
    "GCD Sum",
    "DFT",
    "1D Convolution",
    "Dot Product",
    "Harmonic Sum",
    "Bubble Sort",
    "Kadane's Algorithm",
  ];

  // Array of WASM functions (imported)
  const wasmFunctions = [
    algo1, algo2, algo3, algo4, algo5,
    algo6, algo7, algo8, algo9, algo10,
    algo11, algo12, algo13, algo14, algo15,
  ];

  // --- JavaScript Implementations ---

  const algo1_js = () => {
    // Matrix Multiplication
    const size = 400;
    const modVal = 10007;
    const matrixA = new Uint32Array(size * size);
    const matrixB = new Uint32Array(size * size);
    const result = new Uint32Array(size * size);
    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        matrixA[i * size + j] = (i + j) % modVal;
        matrixB[i * size + j] = ((i * j) % modVal);
      }
    }
    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size; j++) {
        let sum = 0;
        for (let k = 0; k < size; k++) {
          sum = (sum + matrixA[i * size + k] * matrixB[k * size + j]) % modVal;
        }
        result[i * size + j] = sum;
      }
    }
    let finalResult = 0;
    for (let i = 0; i < result.length; i++) {
      finalResult = (finalResult + result[i]) % modVal;
    }
    return finalResult;
  };

  const algo2_js = () => {
    const iterations = 10_000_000;
    const modVal = 1_000_000_007;
    let a = 0, b = 1;
    for (let i = 0; i < iterations; i++) {
      const c = (a + b) % modVal;
      a = b;
      b = c;
    }
    return a;
  };

  const algo3_js = () => {
    const limit = 50000;
    const modVal = 10007;
    const isPrime = new Array(limit + 1).fill(true);
    if(limit >= 0) isPrime[0] = false;
    if(limit >= 1) isPrime[1] = false;
    const r = Math.floor(Math.sqrt(limit));
    for (let i = 2; i <= r; i++) {
      if (isPrime[i]) {
        for (let j = i * i; j <= limit; j += i) {
          isPrime[j] = false;
        }
      }
    }
    let sum = 0;
    for (let i = 2; i <= limit; i++) {
      if (isPrime[i]) {
        sum = (sum + i) % modVal;
      }
    }
    return sum;
  };

  const algo4_js = () => {
    const size = 100000;
    const modVal = 10007;
    const arr = new Uint32Array(size);
    let seed = 42;
    for (let i = 0; i < size; i++) {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      arr[i] = seed % 100000;
    }
    const arrSorted = Array.from(arr);
    arrSorted.sort((a, b) => a - b);
    let sum = 0;
    for (let i = 0; i < arrSorted.length; i++) {
      sum = (sum + arrSorted[i]) % modVal;
    }
    return sum;
  };

  const algo5_js = () => {
    const width = 300;
    const height = 300;
    const maxIter = 100;
    let count = 0;
    for (let i = 0; i < width; i++) {
      for (let j = 0; j < height; j++) {
        const cx = (i / width) * 3.5 - 2.5;
        const cy = (j / height) * 2.0 - 1.0;
        let x = 0, y = 0;
        let iter = 0;
        while (x*x + y*y <= 4 && iter < maxIter) {
          const temp = x*x - y*y + cx;
          y = 2*x*y + cy;
          x = temp;
          iter++;
        }
        if (iter === maxIter) {
          count++;
        }
      }
    }
    return count;
  };

  const algo6_js = () => {
    const iterations = 1_000_000;
    let inside = 0;
    let seed = 42;
    for (let i = 0; i < iterations; i++) {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      const x = (seed % 10000) / 10000;
      seed = (seed * 1664525 + 1013904223) >>> 0;
      const y = (seed % 10000) / 10000;
      if (x*x + y*y <= 1) {
        inside++;
      }
    }
    const pi_est = 4 * inside / iterations;
    return Math.floor(pi_est * 1000);
  };

  const algo7_js = () => {
    const n = 10;
    let count = 0;
    const cols = new Array(n).fill(false);
    const diag1 = new Array(2 * n - 1).fill(false);
    const diag2 = new Array(2 * n - 1).fill(false);
    function solve(row: number) {
      if (row === n) {
        count++;
        return;
      }
      for (let col = 0; col < n; col++) {
        if (!cols[col] && !diag1[row + col] && !diag2[row + n - col - 1]) {
          cols[col] = diag1[row + col] = diag2[row + n - col - 1] = true;
          solve(row + 1);
          cols[col] = diag1[row + col] = diag2[row + n - col - 1] = false;
        }
      }
    }
    solve(0);
    return count;
  };

  const algo8_js = () => {
    const n = 20000;
    const modVal = 1000000007;
    let fact = 1;
    for (let i = 1; i <= n; i++) {
      fact = (fact * i) % modVal;
    }
    return fact;
  };

  const algo9_js = () => {
    const limit = 100000;
    const modVal = 10007;
    let sum = 0;
    function gcd(a: number, b: number): number {
      while (b !== 0) {
        const temp = a % b;
        a = b;
        b = temp;
      }
      return a;
    }
    for (let i = 1; i <= limit; i++) {
      const g = gcd(i, i + 100);
      sum = (sum + g) % modVal;
    }
    return sum;
  };

  const algo10_js = () => {
    const n = 256;
    const input = [];
    for (let i = 0; i < n; i++) {
      input.push(Math.sin(i));
    }
    let sum = 0;
    for (let k = 0; k < n; k++) {
      let re = 0, im = 0;
      for (let t = 0; t < n; t++) {
        const angle = 2 * Math.PI * t * k / n;
        re += input[t] * Math.cos(angle);
        im -= input[t] * Math.sin(angle);
      }
      sum += Math.sqrt(re*re + im*im);
    }
    return Math.floor(sum);
  };

  const algo11_js = () => {
    const size = 1000;
    const kernelSize = 50;
    const modVal = 10007;
    const data = new Uint32Array(size);
    const kernel = new Uint32Array(kernelSize);
    for (let i = 0; i < size; i++) {
      data[i] = i % 100;
    }
    for (let i = 0; i < kernelSize; i++) {
      kernel[i] = (i + 1) % 10;
    }
    const result = new Uint32Array(size);
    for (let i = 0; i < size; i++) {
      let sum = 0;
      for (let j = 0; j < kernelSize; j++) {
        if (i >= j) {
          sum += data[i - j] * kernel[j];
        }
      }
      result[i] = sum % modVal;
    }
    let finalResult = 0;
    for (let i = 0; i < result.length; i++) {
      finalResult = (finalResult + result[i]) % modVal;
    }
    return finalResult;
  };

  const algo12_js = () => {
    const size = 1_000_000;
    const modVal = 10007;
    const vec1 = new Uint32Array(size);
    const vec2 = new Uint32Array(size);
    for (let i = 0; i < size; i++) {
      vec1[i] = i % 1000;
      vec2[i] = (size - i) % 1000;
    }
    let dot = 0;
    for (let i = 0; i < size; i++) {
      dot = (dot + vec1[i] * vec2[i]) % modVal;
    }
    return dot;
  };

  const algo13_js = () => {
    const n = 10_000_000;
    let sum = 0;
    for (let i = 1; i <= n; i++) {
      sum += 1 / i;
    }
    return Math.floor(sum * 1000);
  };

  const algo14_js = () => {
    const size = 1000;
    const modVal = 10007;
    const arr = [];
    for (let i = 0; i < size; i++) {
      arr.push(size - i);
    }
    for (let i = 0; i < size; i++) {
      for (let j = 0; j < size - 1; j++) {
        if (arr[j] > arr[j+1]) {
          [arr[j], arr[j+1]] = [arr[j+1], arr[j]];
        }
      }
    }
    let sum = 0;
    for (let i = 0; i < arr.length; i++) {
      sum = (sum + arr[i]) % modVal;
    }
    return sum;
  };

  const algo15_js = () => {
    const size = 100000;
    const arr = new Array<number>(size);
    let seed = 42;
    for (let i = 0; i < size; i++) {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      arr[i] = (seed % 1001) - 500;
    }
    let maxCurrent = arr[0];
    let maxGlobal = arr[0];
    for (let i = 1; i < size; i++) {
      maxCurrent = Math.max(arr[i], maxCurrent + arr[i]);
      if (maxCurrent > maxGlobal) {
        maxGlobal = maxCurrent;
      }
    }
    return maxGlobal < 0 ? 0 : maxGlobal;
  };

  // Array of JS algorithm functions.
  const jsFunctions = [
    algo1_js, algo2_js, algo3_js, algo4_js, algo5_js,
    algo6_js, algo7_js, algo8_js, algo9_js, algo10_js,
    algo11_js, algo12_js, algo13_js, algo14_js, algo15_js,
  ];

  // Run both WASM and JS versions for the given algorithm (by index)
  const runAlgorithm = async (index: number) => {
    await init(); // Ensure WASM is initialized
    const wasmFunc = wasmFunctions[index];
    const startWasm = performance.now();
    const wasmRes = wasmFunc();
    const endWasm = performance.now();
    const jsFunc = jsFunctions[index];
    const startJS = performance.now();
    const jsRes = jsFunc();
    const endJS = performance.now();
    const newResults = [...results];
    newResults[index] = {
      wasmTime: endWasm - startWasm,
      wasmResult: wasmRes,
      jsTime: endJS - startJS,
      jsResult: jsRes
    };
    setResults(newResults);
  };

  return (
    <div style={{ padding: '2rem' }}>
      <h1>Rust vs JS – 15 Algorithms</h1>
      {algoLabels.map((label, index) => (
        <div key={index} style={{ marginBottom: '1rem', border: '1px solid #ccc', padding: '1rem' }}>
          <h2>{index + 1}. {label}</h2>
          <button onClick={() => runAlgorithm(index)}>Run Algorithm {index + 1}</button>
          <div>
            {results[index].wasmTime !== null && (
              <p>
                WASM time: {results[index].wasmTime.toFixed(2)}ms, Result: {results[index].wasmResult}
              </p>
            )}
            {results[index].jsTime !== null && (
              <p>
                JS time: {results[index].jsTime.toFixed(2)}ms, Result: {results[index].jsResult}
              </p>
            )}
          </div>
        </div>
      ))}
    </div>
  );
}

export default App;
