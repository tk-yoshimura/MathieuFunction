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

                    bool is_nan = false;
            
                    for (double q = 0; q <= 8 * Math.Max(1, n * n); q += Math.Max(1, n * n) / 64d) {
                        double x = EigenFP64.InitialValue(func, n, q);
                        double y = EigenFP64.SearchFit(func, n, q, x).value;

                        if (double.IsNaN(y)) {
                            is_nan = true;
                        }

                        sw.WriteLine($"{q},{x},{y}");
                    }

                    if (is_nan) {
                        Console.WriteLine("detected nan");
                    }
                }
            }

            //EigenFunc func = EigenFunc.B;
            //int n = 256;
            //double q = 36864;
            //
            //double x = EigenFP64.InitialValue(func, n, q);
            //double y = EigenFP64.SearchFit(func, n, q, x).value;
            
            Console.WriteLine("END");
            Console.Read();
        }
    }
}