using System;
using System.Xml.Serialization;

namespace MathieuMP {
    class Program {
        static void Main() {
            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
            
                for (int n = func == EigenFunc.A ? 0: 1; n <= 256; n++) {
                    Console.WriteLine($"{func}{n}");

                    using StreamWriter sw = new($"../../../../results/eigen_{func}_{n}_approx.csv");
                    
                    sw.WriteLine("q,approx,convergence");

                    for (double q = 0; q <= 8 * Math.Max(1, n * n); q += Math.Max(1, n * n) / 64d) {
                        double x = EigenFP64.InitialValue(func, n, q);
                        double y = EigenFP64.SearchFit(func, n, q, x).value;
            
                        sw.WriteLine($"{q},{x},{y}");
                    }
                }
            }
            
            Console.WriteLine("END");
            Console.Read();
        }
    }
}