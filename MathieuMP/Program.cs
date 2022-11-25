using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            //for (double q = 0; q <= 64; q += 1d / 64) {                
            //    double a = EigenFP64.InitialValue(EigenFunc.B, 4, q);
            //    double b = EigenFP64.Value(EigenFunc.B, 4, q, zero_shift: true).value;
            //    Console.WriteLine($"{q},{a},{b}");
            //}

            //EigenFP64.Value(EigenFunc.B, 4, 37.375, zero_shift: true);

            //EigenFP64.Fraction(EigenFunc.B, 4, 37.375, -12);

            EigenFP64.Fraction(EigenFunc.B, 4, 37.3753585154754, -12);
            
            for (double a = -11.875; a >= -12.125; a -= 1d / 4096) {
                double d = EigenFP64.Fraction(EigenFunc.B, 4, 37.40, a);
                Console.WriteLine($"{a},{d}");
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}