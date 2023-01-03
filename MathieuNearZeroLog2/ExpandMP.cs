using MultiPrecision;

namespace MathieuNearZero {

    internal struct Plus4<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 4);
    }

    internal struct Plus8<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 8);
    }

    internal struct Plus16<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 16);
    }

    internal struct Plus32<N> : IConstant where N : struct, IConstant {
        public int Value => checked(default(N).Value + 32);
    }
}
