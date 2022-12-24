using MultiPrecision;
using System;

namespace MathieuMP {
    class Program {
        static void Main() {
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                for (int n = func == EigenFunc.A ? 0 : 1; n <= 256; n++) {
                    Console.WriteLine($"{func}{n}");
            
                    using StreamWriter sw = new($"../../../../results/needs_frac_mp4_{func}_{n}.csv");
            
                    int terms = 1;
            
                    for (MultiPrecision<Pow2.N4> q = Math.Max(1, n * n); q <= 64 * Math.Max(1, n * n); q += Math.Max(1, n * n)) {
                        (MultiPrecision<Pow2.N4> a, terms) = EigenMP<Pow2.N4>.ConvergenceFracTerms(func, n, q, Math.Max(1, terms - 2));
                        
                        Console.WriteLine($"{q},{a},{terms}");
                        sw.WriteLine($"{q},{a},{terms}");
                    }
                }
            }
            
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
                for (int n = func == EigenFunc.A ? 0 : 1; n <= 256; n++) {
                    Console.WriteLine($"{func}{n}");
            
                    using StreamWriter sw = new($"../../../../results/needs_frac_mp8_{func}_{n}.csv");
            
                    int terms = 1;
            
                    for (MultiPrecision<Pow2.N8> q = Math.Max(1, n * n); q <= 64 * Math.Max(1, n * n); q += Math.Max(1, n * n)) {
                        (MultiPrecision<Pow2.N8> a, terms) = EigenMP<Pow2.N8>.ConvergenceFracTerms(func, n, q, Math.Max(1, terms - 2));
                        
                        //Console.WriteLine($"{q},{a},{terms}");
                        sw.WriteLine($"{q},{a},{terms}");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}