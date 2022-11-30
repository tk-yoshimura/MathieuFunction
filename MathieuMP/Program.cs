using System.Xml.Serialization;

namespace MathieuMP {
    class Program {
        static void Main() {
            //for (double q = 0; q <= 64; q += 1d / 64) {
            //    double a = EigenFP64.InitialValue(EigenFunc.B, 8, q);
            //    double b = EigenFP64.Value(EigenFunc.B, 8, q, zero_shift: true).value;
            //    Console.WriteLine($"{q},{a},{b}");
            //}

            //static double f(double x) {
            //    return EigenFP64.Fraction(EigenFunc.B, 4, 37.3754, x);
            //}
            //
            //for (double x = -24; x <= 16; x += 1d / 32) {
            //    double y = f(x);
            //
            //    bool b = SequenceUtil.IsReciprocalCurve(f(x - 3d / 64), f(x - 2d / 64), f(x - 1d / 64), f(x), f(x + 1d / 64), f(x + 2d / 64), f(x + 3d / 64));
            //
            //    Console.WriteLine($"{x},{y},{(b ? 1 : 0)}");
            //}
            
            EigenFP64.SearchFit(EigenFunc.B, 4, 37.375, -12);
            EigenFP64.SearchFit(EigenFunc.B, 4, 37.3754, -12.001);
            EigenFP64.SearchFit(EigenFunc.B, 4, 37.3753585154754, -12);

            //double a = EigenFP64.Fraction(EigenFunc.B, 4, 37.3754, -12.00003504803275);

           //for (double a = -24; a <= 16; a += 1d / 32) {
           //    double d = EigenFP64.Fraction(EigenFunc.B, 8, 37.25, a, sqrt_scaling: true);
           //    Console.WriteLine($"{a},{d}");
           //}

            //bool islinear = EigenFP64.IsLinear(0, 1, 2, 3.8, 4);

            //for (double p = 0.1; p <= 2.0; p += 0.1) {
            //    (double x, double error) = RootFinder.SecantSearch((x) => Math.Cos(x) - 0.4, p, 0.1);
            //
            //    Console.WriteLine($"{x}, {error}\n");
            //}
            //
            //for (double p = 0.1; p <= 2.0; p += 0.1) {
            //    (double x, double error) = RootFinder.SecantSearch((x) => Math.Sin(x) - 0.4, p, 0.1);
            //
            //    Console.WriteLine($"{x}, {error}\n");
            //}

            Console.WriteLine("END");
            Console.Read();
        }
    }
}