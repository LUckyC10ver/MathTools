using System;
using System.IO;
using MathNet.Numerics.LinearAlgebra;

namespace MathTools.Core.Legacy.LinearAlgebra
{
    /// <summary>
    /// ILUT 预条件器的简化实现（基于 LU 分解）。
    /// </summary>
    public sealed class BCILUT
    {
        private Matrix<double> _lower;
        private Matrix<double> _upper;
        private int _dropped;
        private double _eps;
        private int _k;

        public BCILUT(Matrix<double> matrix, int k, double eps)
        {
            if (matrix == null)
            {
                throw new ArgumentNullException(nameof(matrix));
            }

            _k = k;
            _eps = eps;
            Init(matrix);
        }

        public int Dropped => _dropped;
        public Matrix<double> Lower => _lower;
        public Matrix<double> Upper => _upper;

        public void Init(Matrix<double> matrix)
        {
            var lu = matrix.LU();
            _lower = lu.L;
            _upper = lu.U;
            _dropped = 0;
        }

        public Vector<double> Solve(Vector<double> rhs)
        {
            if (_lower == null || _upper == null)
            {
                throw new InvalidOperationException("ILUT not initialized");
            }

            var y = _lower.Solve(rhs);
            return _upper.Solve(y);
        }

        public void Output(TextWriter writer)
        {
            if (writer == null)
            {
                throw new ArgumentNullException(nameof(writer));
            }

            writer.WriteLine("BCILUT");
            writer.WriteLine($"Dropped: {_dropped}");
            writer.WriteLine($"Threshold: {_eps}");
            writer.WriteLine($"Fill depth: {_k}");
        }
    }
}
