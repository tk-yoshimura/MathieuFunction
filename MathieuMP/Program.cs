using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            for (double q = 0; q <= 64; q += 1d / 64) {                
                double a = EigenFP64.InitialValue(EigenFunc.B, 4, q);
                double b = EigenFP64.Value(EigenFunc.B, 4, q, zero_shift: true).value;
                Console.WriteLine($"{q},{a},{b}");
            }

            //for (double q = 0; q <= 64; q += 1d / 64) {              
            //    double a = EigenFP64.InitialValue(EigenFunc.B, 3, q);
            //    Console.WriteLine($"{q},{a}");
            //}

            Console.WriteLine("END");
            Console.Read();
        }
    }
}