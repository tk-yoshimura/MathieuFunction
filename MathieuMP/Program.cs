using System;
using System.Xml.Serialization;

namespace MathieuMP {
    class Program {
        static void Main() {

            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.375358, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.375359, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.37535851, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.37535852, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.37535851547, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.3753585154754, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.37535851548, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.375, -11.05);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.375, -12.05);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.375, -12);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.3754, -12.001);

            foreach (EigenFunc func in new[] { EigenFunc.A, EigenFunc.B }) {
            
                for (int n = func == EigenFunc.A ? 0: 1; n <= 256; n++) {
                    Console.WriteLine($"{func}{n}");
            
                    using StreamWriter sw = new($"../../../../results/eigen_{func}_{n}_approx.csv");
                    
                    sw.WriteLine("q,approx,convergence,score");
            
                    bool is_nan = false;
            
                    for (double q = 0; q <= 8 * Math.Max(1, n * n); q += Math.Max(1, n * n) / 4096d) {
                        double x = EigenFP64.InitialValue(func, n, q);
                        (double y, double score, _) = EigenFP64.SearchFit(func, n, q, x);
            
                        if (double.IsNaN(y)) {
                            is_nan = true;
                        }
            
                        sw.WriteLine($"{q},{x},{y},{score}");
                    }
            
                    if (is_nan) {
                        Console.WriteLine("detected nan");
                    }
                }
            }

            //EigenFunc func = EigenFunc.B;
            //int n = 162;
            //double q = 20791.4501953125;
            //
            //double x = EigenFP64.InitialValue(func, n, q);
            //double y = EigenFP64.SearchFit(func, n, q, x).value;
            //
            //for (double a = 9000; a <= 10000; a += Math.ScaleB(1, -2)) {
            //    double d = EigenFP64.Fraction(func, n, q, a);
            //    Console.WriteLine($"{a},{d}");
            //}
            //
            //double score6 = RootFinder.LinearityScore((x) => 1 / x, 0);
            //double score1 = RootFinder.LinearityScore((x) => 2 * x, 0);
            //double score2 = RootFinder.LinearityScore((x) => x, 0);
            //double score3 = RootFinder.LinearityScore((x) => x - 2, 2);
            //double score4 = RootFinder.LinearityScore((x) => x / 256, 0);
            //double score5 = RootFinder.LinearityScore((x) => x * x, 0);
            //double score7 = RootFinder.LinearityScore((x) => x - 2 + 1e-12, 2);
            //double score8 = RootFinder.LinearityScore((x) => 1 / (x - 2) + 1e-12, 2);
            //double score9 = RootFinder.LinearityScore((x) => 1 / (x - 2) + 1e-12, 1.99999999);
            //double score10 = RootFinder.LinearityScore((x) => 1 / (x - 2) + 1e-12, 1.99999999999999);
            
            Console.WriteLine("END");
            Console.Read();
        }
    }
}