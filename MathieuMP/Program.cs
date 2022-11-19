using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            for (double a = -16; a <= 16; a += 1d / 32) {
                double d = MathieuEigenFP64.FractionA(3, 4, a);

                Console.WriteLine($"{a},{d}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}