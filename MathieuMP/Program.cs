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

            EigenFP64.SearchFit(EigenFunc.B, 4, 37.3754, -12.001);
            //EigenFP64.SearchFit(EigenFunc.B, 4, 37.3753585154754, -12);

            //double a = EigenFP64.Fraction(EigenFunc.B, 4, 37.375, -11.9996337890625);

            //for (double a = -12 + 1d/64; a >= -12 - 1d/64; a -= 1d / 65536) {
            //    double d = EigenFP64.Fraction(EigenFunc.B, 4, 37.3754, a);
            //    Console.WriteLine($"{a},{d}");
            //}

            //bool islinear = EigenFP64.IsLinear(0, 1, 2, 3.8, 4);

            Console.WriteLine("END");
            Console.Read();
        }
    }
}