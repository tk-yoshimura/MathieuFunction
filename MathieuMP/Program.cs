using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            //for (double a = -16; a <= 16; a += 1d / 32) {
            //    double d = MathieuEigenFP64.FractionB(3, 4, a);
            //
            //    Console.WriteLine($"{a},{d}");
            //}

            for (int n = 1; n <= 16; n++) {
                (double y, double err) = MathieuEigenFP64.SearchFitB(n, Math.ScaleB(1, -400), 0);
            
                Console.WriteLine(y);
                Console.WriteLine(err);
            }


            Console.WriteLine("END");
            Console.Read();
        }
    }
}