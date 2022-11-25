using MultiPrecision;

namespace MathieuMP {
    class Program {
        static void Main() {
            for (int a = -2; a <= 2; a++) { 
                for (int b = -2; b <= 2; b++) {
                    for (int c = -2; c <= 2; c++) {
                        double x = EigenFP64.RootCubic(a, b, c);
                        double z = a + b * x + c * x * x + x * x * x;
                        Console.WriteLine(x);
                        Console.WriteLine(z);       
                        Console.WriteLine("");
                    }
                }
            }

            Console.WriteLine("END");
            Console.Read();
        }
    }
}