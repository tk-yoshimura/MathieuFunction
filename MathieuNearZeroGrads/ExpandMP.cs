using MultiPrecision;

namespace MathieuMP {

    internal struct N80 : IConstant {
        public int Value => 80;
    }

    internal struct Plus4<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 4);
    }

    internal struct Plus8<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 8);
    }

    internal struct Plus16<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 16);
    }

    internal struct Plus24<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 24);
    }

    internal struct Plus32<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 32);
    }

    internal struct Plus40<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 40);
    }

    internal struct Plus48<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 48);
    }

    internal struct Plus56<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 56);
    }
}
