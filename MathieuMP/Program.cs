using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {

            //for (int n = 0; n <= 16; n++) {
            //    Console.WriteLine($"n={n}");
            //    for (double q = 0, h = 16 * Math.Max(1, n * n); q <= h; q += h / 64) {
            //        double init = EigenFP64.InitialValue(EigenFunc.A, n, q);
            //        (double value, double err) = EigenFP64.Value(EigenFunc.A, n, q, zero_shift: true);
            //
            //        Console.WriteLine($"{q},{init},{value},{err}");
            //    }
            //}

            EigenFP64.Value(EigenFunc.A, 4, 15, zero_shift: true);
            Console.WriteLine("");

            EigenFP64.Value(EigenFunc.A, 4, 16, zero_shift: true);
            Console.WriteLine("");

            EigenFP64.Value(EigenFunc.A, 4, 17, zero_shift: true);
            Console.WriteLine("");

            EigenFP64.Value(EigenFunc.A, 4, 18, zero_shift: true);
            Console.WriteLine("");

            for (double a = 0; a <= 16; a += 1d / 64) {
                double y = EigenFP64.Fraction(EigenFunc.A, 4, 16, a);
                Console.WriteLine($"{a},{y}");
            }


            Console.WriteLine("END");
            Console.Read();
        }
    }
}